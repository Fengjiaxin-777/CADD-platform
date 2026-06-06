# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
import io
import sys

# 解决输出流编码问题
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

st.set_page_config(page_title="数据质控与清洗", layout="wide")

# 统一注入样式
st.markdown("""
<style>
    .main { background-color: #f8fafc; }
    .blue-banner {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
        color: #ffffff;
        padding: 30px;
        border-radius: 12px;
        margin-bottom: 25px;
    }
    .blue-banner h1 { color: #ffffff !important; margin: 0 0 8px 0 !important; font-size: 28px; }
    .blue-banner p { color: #cbd5e1 !important; margin: 0 !important; font-size: 14px; }
    .card {
        background-color: #ffffff;
        padding: 24px;
        border-radius: 12px;
        border: 1px solid #e2e8f0;
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05);
        margin-bottom: 20px;
    }
    .section-header {
        background-color: #f1f5f9;
        padding: 10px 16px;
        border-radius: 8px;
        color: #0f172a;
        font-weight: 700;
        font-size: 16px;
        margin-bottom: 18px;
        border-left: 5px solid #1e3a8a;
        display: flex;
        align-items: center;
    }
    .report-card {
        background-color: #ffffff;
        border-left: 4px solid #1e3a8a;
        padding: 20px;
        border-radius: 8px;
        border: 1px solid #e2e8f0;
        margin-top: 15px;
    }
    .stButton>button {
        background-color: #1e3a8a !important;
        color: white !important;
        border-radius: 8px !important;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>分子数据质控与清洗中心</h1>
    <p>提供自定义阈值的理化性质筛选通道。清洗后的数据将直接作为缓存，供应给下游的物性剖析、药代评估和模型训练模块。</p>
</div>
""", unsafe_allow_html=True)

# 初始化计算所得的数据容器
if "raw_calculated_df" not in st.session_state:
    st.session_state["raw_calculated_df"] = None

# ==========================================
# 布局配置 1：文件上传与格式向导
# ==========================================
st.markdown('<div class="card">', unsafe_allow_html=True)
st.markdown('<div class="section-header">导入分子结构原始文件</div>', unsafe_allow_html=True)

col_up, col_tips = st.columns([0.45, 0.55])

with col_up:
    uploaded_file = st.file_uploader("选择待质控的 CSV 表格文件：", type=["csv"])
    
with col_tips:
    st.markdown("""
    <div style="background-color: #f8fafc; border: 1px dashed #cbd5e1; padding: 15px; border-radius: 8px; font-size: 13px; color: #475569; line-height: 1.6;">
        <b>⚠️ 原始文件表头规范说明：</b><br>
        1. <b>核心主键</b>：必须含有名为 <code>smiles</code> 的列，系统将其作为判定分子合法性及三维拓扑特征的基础。<br>
        2. <b>其他特征（选填）</b>：可以包含 <code>label</code> 或其他活性测试字段，清洗时这些数据会同步保留输出。
    </div>
    """, unsafe_allow_html=True)

st.markdown('</div>', unsafe_allow_html=True)

# ==========================================
# 理化性质全局计算（只在初次载入或重换文件时运算）
# ==========================================
if uploaded_file is not None:
    # 读入原始表格
    df_input = pd.read_csv(uploaded_file)
    
    if "smiles" not in df_input.columns:
        st.error("表格结构异常：未能在您的数据表中检测到 'smiles' 字段列，请检查表头名称。")
    else:
        # 检测是否需要为当前文件重新计算描述符
        # 如果缓存为空，或者上传的文件行数/内容与缓存内容不合，则触发计算
        trigger_calc = False
        if st.session_state["raw_calculated_df"] is None:
            trigger_calc = True
        else:
            if len(st.session_state["raw_calculated_df"]) != len(df_input):
                trigger_calc = True
                
        if trigger_calc:
            with st.spinner("正在解析化学结构并预计算分子描述符参数..."):
                mws, logps, tpsas, hbds, hbas, rbs, valids = [], [], [], [], [], [], []
                
                for s in df_input["smiles"]:
                    try:
                        s_str = str(s).strip()
                        m = Chem.MolFromSmiles(s_str)
                        if m:
                            mws.append(Descriptors.MolWt(m))
                            logps.append(Crippen.MolLogP(m))
                            tpsas.append(Descriptors.TPSA(m))
                            hbds.append(Lipinski.NumHDonors(m))
                            hbas.append(Lipinski.NumHAcceptors(m))
                            rbs.append(Descriptors.NumRotatableBonds(m))
                            valids.append(True)
                        else:
                            # 异常 SMILES
                            mws.append(np.nan); logps.append(np.nan); tpsas.append(np.nan)
                            hbds.append(np.nan); hbas.append(np.nan); rbs.append(np.nan)
                            valids.append(False)
                    except:
                        mws.append(np.nan); logps.append(np.nan); tpsas.append(np.nan)
                        hbds.append(np.nan); hbas.append(np.nan); rbs.append(np.nan)
                        valids.append(False)
                        
                df_calc = df_input.copy()
                df_calc["MW"] = mws
                df_calc["LogP"] = logps
                df_calc["TPSA"] = tpsas
                df_calc["HBD"] = hbds
                df_calc["HBA"] = hbas
                df_calc["RotatableBonds"] = rbs
                df_calc["_is_valid_smiles_"] = valids
                
                st.session_state["raw_calculated_df"] = df_calc
                st.success("理化参数初始化转换完毕。")

# ==========================================
# 布局配置 2：理化限制区间自主配置（滑块过滤）
# ==========================================
if st.session_state["raw_calculated_df"] is not None:
    df_work = st.session_state["raw_calculated_df"].copy()
    
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown('<div class="section-header">自定义分子性质质控阈值范围</div>', unsafe_allow_html=True)
    
    st.markdown("*调节以下各物理化学常数的允许区间，系统将依据您选出的规则过滤出高纯净度的数据载体。*")
    
    col_filter_1, col_filter_2, col_filter_3 = st.columns(3)
    
    with col_filter_1:
        # 分子量过滤区间
        mw_min, mw_max = float(df_work["MW"].min()), float(df_work["MW"].max())
        if np.isnan(mw_min) or np.isnan(mw_max) or mw_min == mw_max:
            mw_min, mw_max = 0.0, 1000.0
        mw_range = st.slider(
            "分子量允许范围 (MW, Da)：",
            min_value=0.0,
            max_value=1200.0,
            value=(max(0.0, mw_min), min(800.0, mw_max)),
            step=10.0
        )
        
        # 脂水分配系数过滤区间
        logp_min, logp_max = float(df_work["LogP"].min()), float(df_work["LogP"].max())
        if np.isnan(logp_min) or np.isnan(logp_max) or logp_min == logp_max:
            logp_min, logp_max = -5.0, 10.0
        logp_range = st.slider(
            "脂水分配系数允许范围 (LogP)：",
            min_value=-8.0,
            max_value=12.0,
            value=(max(-4.0, logp_min), min(6.0, logp_max)),
            step=0.1
        )
        
    with col_filter_2:
        # 极性拓扑表面积过滤区间
        tpsa_min, tpsa_max = float(df_work["TPSA"].min()), float(df_work["TPSA"].max())
        if np.isnan(tpsa_min) or np.isnan(tpsa_max) or tpsa_min == tpsa_max:
            tpsa_min, tpsa_max = 0.0, 300.0
        tpsa_range = st.slider(
            "极性拓扑表面积范围 (TPSA, Å²)：",
            min_value=0.0,
            max_value=350.0,
            value=(max(0.0, tpsa_min), min(150.0, tpsa_max)),
            step=5.0
        )
        
        # 柔性旋转键数过滤区间
        rb_min, rb_max = float(df_work["RotatableBonds"].min()), float(df_work["RotatableBonds"].max())
        if np.isnan(rb_min) or np.isnan(rb_max) or rb_min == rb_max:
            rb_min, rb_max = 0.0, 20.0
        rb_range = st.slider(
            "可旋转单键个数上限和下限 (RotatableBonds)：",
            min_value=0,
            max_value=30,
            value=(int(max(0, rb_min)), int(min(12, rb_max))),
            step=1
        )
        
    with col_filter_3:
        # 氢键供体数范围
        hbd_min, hbd_max = float(df_work["HBD"].min()), float(df_work["HBD"].max())
        if np.isnan(hbd_min) or np.isnan(hbd_max) or hbd_min == hbd_max:
            hbd_min, hbd_max = 0.0, 15.0
        hbd_range = st.slider(
            "氢键供体数量范围 (HBD)：",
            min_value=0,
            max_value=20,
            value=(int(max(0, hbd_min)), int(min(6, hbd_max))),
            step=1
        )
        
        # 氢键受体数范围
        hba_min, hba_max = float(df_work["HBA"].min()), float(df_work["HBA"].max())
        if np.isnan(hba_min) or np.isnan(hba_max) or hba_min == hba_max:
            hba_min, hba_max = 0.0, 25.0
        hba_range = st.slider(
            "氢键受体数量范围 (HBA)：",
            min_value=0,
            max_value=25,
            value=(int(max(0, hba_min)), int(min(12, hba_max))),
            step=1
        )
        
    # 其他常用清洗项选项
    st.markdown("---")
    st.markdown("**结构规范性与除杂处理选项：**")
    col_opt1, col_opt2 = st.columns(2)
    with col_opt1:
        drop_duplicates = st.checkbox("去重：自动剔除包含重复 SMILES 的化合物记录", value=True)
    with col_opt2:
        drop_invalid = st.checkbox("质量审计：自动过滤无法被 RDKit 解析的错误化学式结构", value=True)
        
    st.markdown('</div>', unsafe_allow_html=True)
    
    # ==========================================
    # 执行过滤动作并更新全局缓存
    # ==========================================
    # 执行过滤条件
    filtered_df = df_work.copy()
    
    # 1. 踢除无效SMILES
    if drop_invalid:
        filtered_df = filtered_df[filtered_df["_is_valid_smiles_"] == True]
        
    # 2. 去重
    if drop_duplicates:
        filtered_df = filtered_df.drop_duplicates(subset=["smiles"])
        
    # 3. 范围截断
    filtered_df = filtered_df[
        (filtered_df["MW"] >= mw_range[0]) & (filtered_df["MW"] <= mw_range[1]) &
        (filtered_df["LogP"] >= logp_range[0]) & (filtered_df["LogP"] <= logp_range[1]) &
        (filtered_df["TPSA"] >= tpsa_range[0]) & (filtered_df["TPSA"] <= tpsa_range[1]) &
        (filtered_df["RotatableBonds"] >= rb_range[0]) & (filtered_df["RotatableBonds"] <= rb_range[1]) &
        (filtered_df["HBD"] >= hbd_range[0]) & (filtered_df["HBD"] <= hbd_range[1]) &
        (filtered_df["HBA"] >= hba_range[0]) & (filtered_df["HBA"] <= hba_range[1])
    ]
    
    # 清理掉用于标记内部结构的临时字段
    if "_is_valid_smiles_" in filtered_df.columns:
        filtered_df = filtered_df.drop(columns=["_is_valid_smiles_"])
        
    # 保存结果并提交至 Session state
    st.session_state["cleaned_df"] = filtered_df.copy()
    
    # ==========================================
    # 预处理清洗报告明细展示
    # ==========================================
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown('<div class="section-header">数据清洗质控报告</div>', unsafe_allow_html=True)
    
    col_out_l, col_out_r = st.columns([0.4, 0.6])
    
    with col_out_l:
        st.markdown(f"""
        <div class="report-card" style="margin-top:0px; background-color: #f8fafc;">
            <p style="margin: 0; font-size: 14px; color: #1e3a8a; font-weight: bold;">📊 筛选前/后样本统计对比：</p>
            <hr style="border: 0; border-top: 1px solid #cbd5e1; margin: 10px 0;">
            <table style="width:100%; font-size:13px; color:#475569;">
                <tr><td>原始总行数：</td><td style="text-align:right; font-weight:bold;">{len(df_input)} 个</td></tr>
                <tr><td>清洗后保留值：</td><td style="text-align:right; font-weight:bold; color: #16a34a;">{len(filtered_df)} 个</td></tr>
                <tr><td>过滤淘汰行数：</td><td style="text-align:right; font-weight:bold; color: #dc2626;">{len(df_input) - len(filtered_df)} 个</td></tr>
                <tr><td>有效保留率：</td><td style="text-align:right; font-weight:bold;">{(len(filtered_df)/len(df_input)*100):.2f}%</td></tr>
            </table>
        </div>
        """, unsafe_allow_html=True)
        
    with col_out_r:
        st.markdown("🌐 **清洗筛选后的分子库预览：**")
        st.dataframe(filtered_df.head(10), use_container_width=True)
        
        # 提供导出下载
        csv_clean_data = filtered_df.to_csv(index=False).encode('utf-8-sig')
        st.download_button(
            "💾 下载质控后的物理化学性质分子表 (CSV)",
            csv_clean_data,
            "cleaned_properties_dataset.csv",
            "text/csv",
            use_container_width=True
        )
        
    st.markdown('</div>', unsafe_allow_html=True)
else:
    st.info("提示：请先在栏目上方导入待处理的 CSV 数据表以激活数据过滤组件。")

# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf
import sys
import io

# 解决特定环境下可能导致的输出流编码问题
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# 字体配置，多系统兼容方案（特别针对 Linux/Streamlit Cloud 部署环境加入了 WenQuanYi）
def set_matplot_zh_font():
    plt.rcParams['font.sans-serif'] = [
        'SimHei',             # Windows 黑体
        'Microsoft YaHei',    # Windows 微软雅黑
        'Heiti TC',           # macOS 繁体黑体
        'Arial Unicode MS',   # macOS 兼容中文
        'WenQuanYi Micro Hei',# Linux 下的开源中文字体
        'DejaVu Sans',        # 备用英文
        'sans-serif'
    ]
    plt.rcParams['axes.unicode_minus'] = False # 正常显示负号

set_matplot_zh_font()
      
st.set_page_config(page_title="成药性(ADMET)评估", layout="wide")
      
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
    
    /* 统一的区块标题风格（浅灰底色+深蓝左色条） */
    .section-header {
        background-color: #f1f5f9;
        padding: 10px 16px;
        border-radius: 8px;
        color: #0f172a;
        font-weight: 700;
        font-size: 15px;
        margin-top: 5px;
        margin-bottom: 15px;
        border-left: 5px solid #1e3a8a;
    }
    .report-card {
        background-color: #ffffff;
        border-left: 4px solid #1e3a8a;
        padding: 20px;
        border-radius: 8px;
        border: 1px solid #e2e8f0;
        margin-top: 15px;
        margin-bottom: 15px;
    }
    .metric-badge {
        padding: 8px 16px;
        border-radius: 8px;
        color: white;
        font-weight: 700;
        text-align: center;
        font-size: 14px;
        margin-top: 10px;
        margin-bottom: 15px;
    }
    /* 强力重置 Streamlit 原生 border 容器的边框和内边距，统一风格 */
    div[data-testid="stVerticalBlockBorderLine"] {
        background-color: #ffffff !important;
        border: 1px solid #e2e8f0 !important;
        border-radius: 12px !important;
        padding: 22px !important;
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05) !important;
        margin-bottom: 20px !important;
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
    <h1>成药性与 ADMET 药代特征评估</h1>
    <p>评估预测先导分子的口服吸收效能、血脑屏障跨膜情况及潜在宿主器官毒性，辅助化学家快速筛选分子结构。</p>
</div>
""", unsafe_allow_html=True)
      
def calc_advanced_admet(mol):
    res = {}
    res["分子量(MW)"] = round(Descriptors.MolWt(mol), 2)
    res["LogP"] = round(Crippen.MolLogP(mol), 2)
    res["极性表面积(TPSA)"] = round(MolSurf.TPSA(mol), 2)
    res["氢键供体(HBD)"] = Lipinski.NumHDonors(mol)
    res["氢键受体(HBA)"] = Lipinski.NumHAcceptors(mol)
    res["可旋转键(RB)"] = Lipinski.NumRotatableBonds(mol)
    res["芳香环数"] = Lipinski.NumAromaticRings(mol)
    res["QED成药性指数"] = round(Descriptors.qed(mol), 3)
          
    res["Lipinski规则"] = {
        "分子量 <= 500": res["分子量(MW)"] <= 500,
        "LogP <= 5": res["LogP"] <= 5,
        "氢键供体数 <= 5": res["氢键供体(HBD)"] <= 5,
        "氢键受体数 <= 10": res["氢键受体(HBA)"] <= 10
    }
          
    res["Lipinski五规则违规数"] = sum(not v for v in res["Lipinski规则"].values())
    res["Veber规则违规数"] = sum([res["可旋转键(RB)"] > 10, res["极性表面积(TPSA)"] > 140])
    res["Egan规则违规数"] = sum([res["LogP"] < -1 or res["LogP"] > 6, res["极性表面积(TPSA)"] > 130])
      
    res["胃肠道吸收"] = "高吸收" if res["极性表面积(TPSA)"] < 140 and res["LogP"] > 0 else "低吸收"
    res["血脑屏障渗透性"] = "易穿透" if (res["LogP"] > 2 and res["极性表面积(TPSA)"] < 90) else "难穿透"
    res["肝毒性风险"] = "存在高风险" if (res["芳香环数"] >= 4 or res["LogP"] > 4.5) else "较低风险"
      
    score = 0
    score += 25 if res["Lipinski五规则违规数"] == 0 else max(0, 25 - res["Lipinski五规则违规数"] * 8)
    score += 20 if res["Veber规则违规数"] == 0 else 5
    score += 25 * res["QED成药性指数"]
    score += 15 if res["肝毒性风险"] == "较低风险" else 0
    score += 15 if res["胃肠道吸收"] == "高吸收" else 0
    res["综合成药得分"] = round(score, 1)
      
    if res["综合成药得分"] >= 80:
        res["药效评级"] = "优秀候选分子"
    elif res["综合成药得分"] >= 65:
        res["药效评级"] = "常规（可继续优化）"
    elif res["综合成药得分"] >= 45:
        res["药效评级"] = "受限（改良度较低）"
    else:
        res["药效评级"] = "不宜继续开发"
      
    return res
      
def plot_radar(props):
    scores = {
        "分子量符合度": 100 if props["分子量(MW)"] <= 500 else max(10, 100 - (props["分子量(MW)"]-500)/2),
        "脂溶性符合度": 100 if props["LogP"] <= 5 else max(10, 100 - (props["LogP"]-5)*15),
        "极性区结合力": max(10, min(100, 100 - abs(props["极性表面积(TPSA)"]-80)/1.2)),
        "肠胃道吸收": 100 if props["胃肠道吸收"] == "高吸收" else 30,
        "低肝毒安全性": 100 if props["肝毒性风险"] == "较低风险" else 30,
        "QED定量指数": props["QED成药性指数"] * 100
    }
    vals = list(scores.values())
    vals += vals[:1]
    angles = np.linspace(0, 2*np.pi, 6, endpoint=False).tolist()
    angles += angles[:1]
      
    fig, ax = plt.subplots(figsize=(4.5, 4.5), subplot_kw=dict(polar=True))
    ax.plot(angles, vals, "o-", c="#1e3a8a", lw=2, markersize=5)
    ax.fill(angles, vals, alpha=0.2, c="#1e3a8a")
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(list(scores.keys()), fontsize=9, color="#475569")
    ax.set_ylim(0, 100)
    ax.spines['polar'].set_color('#cbd5e1')
    plt.tight_layout()
    return fig

# ==========================================
# 页面排布 1：配置输入区域
# ==========================================
with st.container(border=True):
    st.markdown('<div class="section-header">模式与输入配置</div>', unsafe_allow_html=True)
    mode = st.radio("请选择评测的输入对象形式：", ["评估单个分子 (SMILES)", "批量处理分子常数表 (CSV)"], horizontal=True)

# ==========================================
# 页面排布 2：根据选择的模式，进行不同的流式结果分析
# ==========================================
if mode == "评估单个分子 (SMILES)":
    with st.container(border=True):
        st.markdown('<div class="section-header">单分子输入评估</div>', unsafe_allow_html=True)
        smi = st.text_input("请输入待评估分子的 SMILES：", value="CC(=O)OC1=CC=CC=C1C(=O)O")
        
    if smi:
        smi_clean = smi.strip()
        m = Chem.MolFromSmiles(smi_clean)
        
        if m is None:
            st.error("输入结构解析失败，请确保您输入的 SMILES 符合标准碳氮价键化学标准。")
        else:
            admet = calc_advanced_admet(m)
            fig_rad = plot_radar(admet)
            
            with st.container(border=True):
                col_res_l, col_res_r = st.columns([0.6, 0.4])
                
                with col_res_l:
                    st.markdown('<div class="section-header">物性评估指标状态</div>', unsafe_allow_html=True)
                    
                    # 双标栏展示核心分数与状态等级
                    metric_col1, metric_col2 = st.columns(2)
                    with metric_col1:
                        st.metric("ADMET 综合评估值", f"{admet['综合成药得分']} / 100")
                    with metric_col2:
                        bg_color = {
                            "优秀候选分子": "#1e3a8a",
                            "常规（可继续优化）": "#475569",
                            "受限（改良度较低）": "#94a3b8",
                            "不宜继续开发": "#cbd5e1"
                        }.get(admet["药效评级"], "#475569")
                        
                        st.markdown(f'<div class="metric-badge" style="background-color: {bg_color}; margin-top:20px;">评估状态：{admet["药效评级"]}</div>', unsafe_allow_html=True)
                        
                    st.markdown("**Lipinski 规则通过详情表格：**")
                    rule_check = []
                    for term, passed in admet["Lipinski规则"].items():
                        rule_check.append({"准则指标": term, "判定值": "正常通过" if passed else "❌ 不合规"})
                    st.table(pd.DataFrame(rule_check))
                        
                with col_res_r:
                    st.markdown('<div class="section-header">药效性质指数量化投影</div>', unsafe_allow_html=True)
                    st.pyplot(fig_rad)
                    plt.close(fig_rad)
                    
            # 评估分析优化建议报告
            st.markdown(f"""
            <div class="report-card">
                <div class="section-header" style="background-color: #eff6ff; border-left: 5px solid #2563eb; margin-bottom: 15px;">成药性评估报告与优化建议</div>
                <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                    当前测试分子的成药性评级为 <strong>{admet['药效评级']}</strong>（QED定量得分为 {admet['QED成药性指数']}）。  
                    预测显示小分子在人体内胃肠道的吸收可能为 <strong>{admet['胃肠道吸收']}</strong>（拓扑极性表面积为 {admet['极性表面积(TPSA)']} Å²）。由于当前有 <strong>{admet['Lipinski五规则违规数']}</strong> 项指标不合规：<br><br>
                    - <em>推荐结构改性方向：</em> 如果分子量偏大或脂溶性（LogP）过高，可能导致水溶性较差，制成口服药片时的利用效果不佳。在下一步结构优化时，建议引入极性基团（例如引入微极性的羟基 -OH、氨基 -NH2，或将末端苯环替换为吡啶环等氮杂环），以降低物理分子排斥力并增强溶解性。
                </p>
            </div>
            """, unsafe_allow_html=True)
else:
    # ==========================================
    # 批量处理模式布局重构
    # ==========================================
    with st.container(border=True):
        st.markdown('<div class="section-header">数据源选择</div>', unsafe_allow_html=True)
        data_source = st.radio(
            "请选择批量评估使用的数据来源：",
            ["使用前序页面清洗的缓存数据", "手动上传新的本地 CSV 数据集"],
            horizontal=True
        )
        
        df_batch = None
        if data_source == "使用前序页面清洗的缓存数据":
            if "cleaned_df" in st.session_state:
                df_batch = st.session_state["cleaned_df"]
                st.success("成功关联前序模块清洗后的物理分子缓存数据。")
            else:
                st.error("缓存中不存在数据。请先执行第一页的数据清洗质控，或选择手动上传本地 CSV 文件。")
        else:
            file = st.file_uploader("外部 CSV 数据集导入", type=["csv"])
            if file:
                df_batch = pd.read_csv(file)
          
    if df_batch is not None:
        with st.container(border=True):
            st.markdown('<div class="section-header">批量分子性质运算与报告</div>', unsafe_allow_html=True)
            
            if "smiles" not in df_batch.columns:
                st.error("表格格式异常：表中缺失了必要的 smiles 结构式列。")
            else:
                st.info(f"数据加载完毕：已准备好 {len(df_batch)} 种分子。点击下方按钮开始批量 ADMET 指标估测。")
                
                if "admet_batch_df" not in st.session_state:
                    st.session_state["admet_batch_df"] = None

                if st.button("🚀 开始批量计算", type="primary", use_container_width=True):
                    with st.spinner("系统后台计算中，请稍候..."):
                        results = []
                        for idx, row in df_batch.iterrows():
                            try:
                                s_val = str(row["smiles"]).strip()
                                if pd.isna(row["smiles"]) or not s_val:
                                    continue
                                m_obj = Chem.MolFromSmiles(s_val)
                                if m_obj:
                                    calc_res = calc_advanced_admet(m_obj)
                                    calc_res["smiles"] = s_val
                                    results.append(calc_res)
                            except:
                                continue
                                     
                        if results:
                            st.session_state["admet_batch_df"] = pd.DataFrame(results)
                        else:
                            st.error("未在数据表中解析出任何可用有效的分子结构公式。")
                  
                # 展示计算结果，若已存在，则保持渲染
                if st.session_state["admet_batch_df"] is not None:
                    df_out = st.session_state["admet_batch_df"]
                    
                    st.markdown("---")
                    st.markdown('<div class="section-header">评估完成的数据集预览 (前15条)</div>', unsafe_allow_html=True)
                    st.dataframe(df_out[["smiles", "综合成药得分", "药效评级", "分子量(MW)", "LogP", "极性表面积(TPSA)"]].head(15), use_container_width=True)
                        
                    csv_raw = df_out.to_csv(index=False).encode('utf-8-sig')
                    st.markdown("<br>", unsafe_allow_html=True)
                    st.download_button(
                        "💾 下载本轮 ADMET 批量评估结果报告 (CSV)",
                        csv_raw,
                        "ADMET_batch_report.csv",
                        use_container_width=True
                    )

import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# 字体配置，确保中文不乱码
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="分子性质探索", layout="wide")

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
    .report-card {
        background-color: #ffffff;
        border-left: 4px solid #1e3a8a;
        padding: 20px;
        border-radius: 8px;
        border: 1px solid #e2e8f0;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>分子物理化学性质探索分析</h1>
    <p>可视化查看分子的各项描述符常数在化学空间中的分布状态，诊断性质之间的线性共线性，筛选合理的物理描述特征。</p>
</div>
""", unsafe_allow_html=True)

st.markdown('<div class="card">', unsafe_allow_html=True)
st.subheader("数据源选择")
data_source = st.radio(
    "指定用于物性分析的表格路径：",
    ["调取前序过滤清洗的数据缓存", "从本地上传新的物理常数 CSV 表"]
)

df = None
if data_source == "使用前序清洗缓存数据" or data_source == "调取前序过滤清洗的数据缓存":
    if "cleaned_df" in st.session_state:
        df = st.session_state["cleaned_df"]
        st.success("成功载入前序清洗页面缓存的数据。")
    else:
        st.error("缓存中不存在数据。请先执行数据质控过程，或者在上方切换为手动上传本地文件。")
else:
    f_up = st.file_uploader("导入物理常数数据集", type=["csv"])
    if f_up:
        df = pd.read_csv(f_up)
st.markdown('</div>', unsafe_allow_html=True)

if df is not None:
    ch_map = {
        "MW": "分子量 (MW)",
        "LogP": "脂水分配系数 (LogP)",
        "TPSA": "极性拓扑表面积 (TPSA)",
        "HBD": "氢键供体数 (HBD)",
        "HBA": "氢键受体数 (HBA)",
        "RotatableBonds": "柔性可旋转键数 (RB)"
    }
    
    needed_cols = list(ch_map.keys())
    missing_cols = [c for c in needed_cols if c not in df.columns]
    
    if missing_cols and "smiles" in df.columns:
        with st.spinner("系统检测到物理性质列缺失，正在启动 RDKit 为您提取结构性质..."):
            mws, logps, tpsas, hbds, hbas, rbs = [], [], [], [], [], []
            for s in df["smiles"]:
                try:
                    m = Chem.MolFromSmiles(s)
                    if m:
                        mws.append(Descriptors.MolWt(m))
                        logps.append(Descriptors.MolLogP(m))
                        tpsas.append(Descriptors.TPSA(m))
                        hbds.append(Lipinski.NumHDonors(m))
                        hbas.append(Lipinski.NumHAcceptors(m))
                        rbs.append(Lipinski.NumRotatableBonds(m))
                    else:
                        mws.append(np.nan); logps.append(np.nan); tpsas.append(np.nan)
                        hbds.append(np.nan); hbas.append(np.nan); rbs.append(np.nan)
                except:
                    mws.append(np.nan); logps.append(np.nan); tpsas.append(np.nan)
                    hbds.append(np.nan); hbas.append(np.nan); rbs.append(np.nan)
            df["MW"] = mws
            df["LogP"] = logps
            df["TPSA"] = tpsas
            df["HBD"] = hbds
            df["HBA"] = hbas
            df["RotatableBonds"] = rbs
            df = df.dropna().reset_index(drop=True)
            st.session_state["cleaned_df"] = df

    col_ctrl, col_view = st.columns([0.35, 0.65])
    
    with col_ctrl:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("分析描述符配置")
        selected = st.multiselect("请选取参与考察的自变量描述符：", needed_cols, default=["MW", "LogP", "TPSA"])
        st.markdown('</div>', unsafe_allow_html=True)
        
        if len(selected) >= 2:
            corr_mat = df[selected].corr()
            co_errors = []
            for i in range(len(corr_mat.columns)):
                for j in range(i):
                    val = corr_mat.iloc[i, j]
                    if abs(val) > 0.75:
                        co_errors.append((corr_mat.columns[i], corr_mat.columns[j], val))
            
            diagnostic_msg = ""
            if co_errors:
                diagnostic_msg += "<strong>警报：物理特征间存在高度共线性问题：</strong><br>"
                for f1, f2, val in co_errors:
                    diagnostic_msg += f"- 参数对 <strong>{f1}</strong> 与 <strong>{f2}</strong> 的相关度在 <strong>{val:.2f}</strong> 行列。<br>"
                diagnostic_msg += """
                <br><em>结构优化建议：</em><br>
                由于两项参数携带显著重复的信息，在后续构建定量预测动力学或活性时，容易干扰回归权重的物理拟合，降低模型的精度。建议在上述强相关配对中，手动剔除其中一个性质参数。
                """
            else:
                diagnostic_msg = """
                <strong>系统性质自检通过：</strong><br>
                当前选定物理特征间的 Pearson 相关性绝对值全数保持在 <strong>0.75 以下</strong>。表明各维度包含的信息相互独立，无需担心信息重叠引起回归精度滑坡。您可以放心地将它们一同用于活性定量评估中。
                """
            st.markdown(f"""
            <div class="report-card">
                <h4 style="margin-top:0px; color:#1e3a8a !important;">相关性与特征冗余诊断</h4>
                <p style="color:#475569; font-size:13px; line-height:1.6; margin:0;">
                    {diagnostic_msg}
                </p>
            </div>
            """, unsafe_allow_html=True)

    with col_view:
        if len(selected) >= 2:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("物理属性在先导物库中的分布直方图")
            feature_pairs = [selected[i:i + 2] for i in range(0, len(selected), 2)]
            
            for pair in feature_pairs:
                cols = st.columns(2)
                for idx, feature in enumerate(pair):
                    with cols[idx]:
                        fig = px.histogram(
                            df, x=feature, template="simple_white",
                            color_discrete_sequence=["#1e3a8a"],
                            labels={feature: ch_map[feature]},
                            title=f"{ch_map[feature]} 的分布直方图"
                        )
                        fig.update_layout(height=280, margin=dict(l=15, r=15, t=35, b=15))
                        st.plotly_chart(fig, use_container_width=True)
            
            st.subheader("特征间 Pearson 参数线性相关矩阵图")
            c_heat_l, c_heat_r = st.columns([0.8, 0.2])
            with c_heat_l:
                fig_corr, ax_corr = plt.subplots(figsize=(5.5, 4))
                # 使用中文标注，配合 font_sans-serif 防乱码
                sns.heatmap(corr_mat, annot=True, cmap="Blues", fmt=".2f", ax=ax_corr, square=True,
                            cbar_kws={"shrink": .7}, linewidths=.5)
                # 转化展示为中文名字
                labels_zh = [ch_map[x] for x in selected]
                ax_corr.set_xticklabels(labels_zh, rotation=45, ha="right", fontsize=9)
                ax_corr.set_yticklabels(labels_zh, rotation=0, fontsize=9)
                ax_corr.set_title("物理参数相关相关热图", fontsize=10)
                st.pyplot(fig_corr)
            st.markdown('</div>', unsafe_allow_html=True)
        else:
            st.warning("提示：请在左侧面板中至少勾选两项物理性质参数。")
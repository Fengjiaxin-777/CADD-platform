import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="分子性质探索", layout="wide")

st.markdown("""
<style>
    .main { background-color: #f8fafc; }
    [data-testid="stSidebar"] { background-color: #0f172a !important; }
    [data-testid="stSidebar"] * { color: #cbd5e1 !important; }
    [data-testid="sidebar-nav-container"] { padding-top: 1.5rem !important; }
    [data-testid="sidebar-nav-item"] {
        padding-top: 16px !important;
        padding-bottom: 16px !important;
        margin-top: 12px !important;
        margin-bottom: 12px !important;
        border-radius: 8px !important;
        font-size: 15px !important;
    }
    [data-testid="sidebar-nav-item-active"] { background-color: #1e3a8a !important; font-weight: 700 !important; }
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
    <h1>分子性质空间分布与相关性分析</h1>
    <p>考察物理化学常用常数的多维分布特征，核验特征是否存在过度相关的冗余线性依赖。</p>
</div>
""", unsafe_allow_html=True)

st.markdown('<div class="card">', unsafe_allow_html=True)
st.subheader("数据源选择")
data_source = st.radio("选择进入分析计算的数据路径：", ["读取前序过滤页面缓存结果", "在本地上传全新 CSV 表格"])
df = None
if data_source == "读取前序过滤页面缓存结果":
    if "cleaned_df" in st.session_state:
        df = st.session_state["cleaned_df"]
        st.success("读取成功：已接入过滤清洗后的系统缓存表。")
    else:
        st.error("系统没有检测到清洗数据缓存。请先进行第一步清洗过滤，或者在上方修改为上传本地 CSV。")
else:
    f_up = st.file_uploader("导入自定义分子属性 CSV 文件", type=["csv"])
    if f_up:
        df = pd.read_csv(f_up)
st.markdown('</div>', unsafe_allow_html=True)

if df is not None:
    ch_map = {
        "MW": "分子量 (MW)",
        "LogP": "脂水分配常数 (LogP)",
        "TPSA": "极性拓扑表面积 (TPSA)",
        "HBD": "氢键供体数 (HBD)",
        "HBA": "氢键受体数 (HBA)",
        "RotatableBonds": "可旋转单键数 (RB)"
    }
    
    needed_cols = list(ch_map.keys())
    missing_cols = [c for c in needed_cols if c not in df.columns]
    
    if missing_cols and "smiles" in df.columns:
        with st.spinner("缺失指定物理特征列，正在通过 RDKit 自动提取..."):
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

    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("分析自变量选择")
    selected = st.multiselect("请选取需要考察的物理常数维度：", needed_cols, default=["MW", "LogP", "TPSA"])
    st.markdown('</div>', unsafe_allow_html=True)

    if len(selected) >= 2:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("物理属性分布直方图")
        for i in range(0, len(selected), 2):
            cols = st.columns(2)
            for idx, feature in enumerate(selected[i:i+2]):
                with cols[idx]:
                    fig = px.histogram(
                        df, x=feature, template="simple_white",
                        color_discrete_sequence=["#1e3a8a"],
                        labels={feature: ch_map[feature]},
                        title=f"{ch_map[feature]} 分布频率"
                    )
                    fig.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig, use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("特征自相关分析与评价")
        
        col_heat, col_diag = st.columns([0.5, 0.5])
        corr_mat = df[selected].corr()
        
        with col_heat:
            fig_corr, ax_corr = plt.subplots(figsize=(5, 3.8))
            sns.heatmap(corr_mat, annot=True, cmap="Blues", fmt=".2f", ax=ax_corr, square=True, linewidths=.5)
            labels_zh = [ch_map[x] for x in selected]
            ax_corr.set_xticklabels(labels_zh, rotation=30, ha="right", fontsize=9)
            ax_corr.set_yticklabels(labels_zh, rotation=0, fontsize=9)
            ax_corr.set_title("物理参数 Pearson 相关系数矩阵", fontsize=10)
            st.pyplot(fig_corr)
            
        with col_diag:
            co_errors = []
            for i in range(len(corr_mat.columns)):
                for j in range(i):
                    val = corr_mat.iloc[i, j]
                    if abs(val) > 0.75:
                        co_errors.append((corr_mat.columns[i], corr_mat.columns[j], val))
            
            diagnostic_msg = ""
            if co_errors:
                diagnostic_msg += "<strong>警报：当前自变量分析发现高度强共线性问题：</strong><br>"
                for f1, f2, val in co_errors:
                    diagnostic_msg += f"- 参数对 [<strong>{ch_map[f1]}</strong>] 与 [<strong>{ch_map[f2]}</strong>] 的共线性关联度极高 ({val:.2f})。<br>"
                diagnostic_msg += """
                <br><em>结构优化改进提示：</em><br>
                上述自变量彼此冗余可能会过度干扰回归分析拟合权重，引致泛化能力变差。推荐在实施 QSAR 回归预测建模前，剔除两者中的其中一维度性。
                """
            else:
                diagnostic_msg = """
                <strong>系统互斥分析达标：</strong><br>
                所有选定维度的 Pearson 相关系数全部处于安全水平（绝对值 <strong>&lt; 0.75</strong>）。数据内部多维差异度良好，不会产生多重参数冗余危险。
                """
            st.markdown(f"""
            <div class="report-card" style="height: 100%;">
                <h4 style="margin-top:0px; color:#1e3a8a !important;">共线性与冗余特征诊断结论</h4>
                <p style="color:#475569; font-size:13px; line-height:1.6; margin:0;">
                    {diagnostic_msg}
                </p>
            </div>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    else:
        st.warning("提示：请在上方自变量列表中至少勾选两项物理特征。")

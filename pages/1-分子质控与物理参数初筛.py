import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="数据清洗与质控", layout="wide")

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
        margin-top: 15px;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>数据清洗与质控过滤</h1>
    <p>检验输入化合物库分子式 SMILES 的物理拓扑正确性，剔除分子量或水油常数极度异常的杂质。</p>
</div>
""", unsafe_allow_html=True)

st.markdown('<div class="card">', unsafe_allow_html=True)
st.subheader("控制配置台")
c_up, c_par1, c_par2 = st.columns([0.4, 0.3, 0.3])
with c_up:
    csv_file = st.file_uploader("导入包含 SMILES 的配体数据集 (CSV)", type=["csv"])
with c_par1:
    f_lip = st.checkbox("强制限制 Lipinski 五规则（允许 1 项超标）", value=True)
    f_mw = st.checkbox("限定分子量处于 160 ~ 500 Da 区间", value=True)
with c_par2:
    f_logp = st.checkbox("限定油水分配系数 LogP 处于 -2 ~ 5 之间", value=True)
st.markdown('</div>', unsafe_allow_html=True)

if csv_file:
    df = pd.read_csv(csv_file)
    
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("待质控初始数据预览")
    st.dataframe(df.head(6), use_container_width=True)
    
    btn_start = st.button("启动数据清洗管线", type="primary", use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    if btn_start:
        if "smiles" not in df.columns:
            st.error("数据表字段中缺少名为 'smiles' 的列！")
        else:
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            cleaned_indices = []
            orig_mws, final_mws = [], []
            raw_len = len(df)
            
            for idx, row in df.iterrows():
                if idx % max(1, raw_len // 20) == 0:
                    val = (idx + 1) / raw_len
                    progress_bar.progress(val)
                    status_text.text(f"质控计算中：{int(val*100)}% ({idx+1}/{raw_len})")
                
                s = row["smiles"]
                try:
                    mol = Chem.MolFromSmiles(s)
                    if mol is None: continue
                    
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    orig_mws.append(mw)
                    
                    hbd = Lipinski.NumHDonors(mol)
                    hba = Lipinski.NumHAcceptors(mol)
                    
                    if f_lip:
                        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
                        if violations > 1: continue
                    if f_mw and not (160 <= mw <= 500): continue
                    if f_logp and not (-2 <= logp <= 5): continue
                    
                    cleaned_indices.append(idx)
                    final_mws.append(mw)
                except:
                    continue
            
            progress_bar.empty()
            status_text.empty()
            
            cleaned_df = df.iloc[cleaned_indices].reset_index(drop=True)
            st.session_state["cleaned_df"] = cleaned_df
            
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("质控清洗完毕结果")
            
            col_tbl, col_fig = st.columns([0.6, 0.4])
            with col_tbl:
                st.write("**过滤保留的标准化合物结果 (前10行)：**")
                st.dataframe(cleaned_df.head(10), use_container_width=True)
            with col_fig:
                st.write("**质控前后化合物总数对比图：**")
                fig, ax = plt.subplots(figsize=(4.5, 3.5))
                ax.bar(["质控前原始数", "质控后保留数"], [raw_len, len(cleaned_df)], color=["#cbd5e1", "#1e3a8a"], width=0.45)
                ax.set_ylabel("小分子样本个数")
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                st.pyplot(fig)
            
            st.markdown(f"""
            <div class="report-card">
                <h4 style="margin-top:0px; color:#1e3a8a !important;">清洗过滤诊断报告</h4>
                <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                    系统从小分子配体库的 {raw_len} 个化合物中，筛除掉了 <strong>{raw_len - len(cleaned_df)}</strong> 个不合规分子。
                    最终数据集保留通过率为：<strong>{(len(cleaned_df)/raw_len * 100):.2f}%</strong>。
                    通过质控过滤能够有效规避极端物理化学常数产生的统计过拟合。
                </p>
            </div>
            """, unsafe_allow_html=True)
            
            csv = cleaned_df.to_csv(index=False).encode('utf-8-sig')
            st.write("")
            st.download_button("导出此质控清洗后的 CSV 结果表", csv, "Cleaned_QC_Molecular_Library.csv", use_container_width=True)
            st.markdown('</div>', unsafe_allow_html=True)
else:
    st.info("提示：请先在上方控制台中载入初始 CSV 数据集。")

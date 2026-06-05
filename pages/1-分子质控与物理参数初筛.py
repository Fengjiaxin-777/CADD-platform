import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# 安全学术字体，规避中文乱码报错
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="数据清洗与质控", layout="wide")

st.markdown("""
<style>
    [data-testid="stSidebar"] { background-color: #0f172a !important; }
    [data-testid="stSidebar"] * { color: #e2e8f0 !important; }
    [data-testid="sidebar-nav-container"] { padding-top: 2rem !important; }
    [data-testid="sidebar-nav-item"] { padding-top: 15px !important; padding-bottom: 15px !important; margin: 10px 0 !important; border-radius: 8px !important; }
    [data-testid="sidebar-nav-item-active"] { background-color: #1e3a8a !important; font-weight: 700 !important; }
    .card { background-color: #ffffff; padding: 24px; border-radius: 12px; border: 1px solid #e2e8f0; margin-bottom: 20px; }
    .report-card { background-color: #ffffff; border-left: 4px solid #1e3a8a; padding: 20px; border-radius: 8px; border: 1px solid #e2e8f0; margin-top: 15px; }
</style>
""", unsafe_allow_html=True)

st.header("SMILES 数据清洗与质控过滤")

# 内置高保真测试集（包含典型的药用小分子，利于测试后续 QSAR 与性质多维度特征）
demo_data = {
    "smiles": [
        "CC(=O)OC1=CC=CC=C1C(=O)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "CC1=CN=C(R)N1", 
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C1=CC=C(C=C1)O", "CCO", "CCCCCC", 
        "CC(=O)NC1=CC=C(O)C=C1", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", 
        "CN(C)C1=CC=C(C=C1)C(=O)C2=CC=CC=C2", "CCCC1=NN(C2=C1N=C(NC2=O)CC)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC",
        "CC12C(CC(O)C3C1CCC4(C)C(CC(O)C23)OC(=O)C4)C(=Cc5cc(=O)oc5)C", "C1COCCO1",
        "CCCCCCCCCCCC(=O)O", "CCN(CC)CCNC(=O)C1=CC=C(N)C=C1", "CC(C)(C)NCC(O)C1=CC(O)=C(O)C=C1",
        "CS(=O)(=O)C1=CC=C(CC2=C1N(C(=O)C3=CC=C(Cl)C=C3)C4=C2C=C(OC)C=C4)C",
        "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O", "CC(=O)N", "CI", 
        "C1=CC=CC=C1", "C2=CC=CC=C2", "C3=CC=CC=C3", "C=O", "CN", "CNC",
        "CC(=O)O", "CCOCC", "CCCCCO", "C1CCCCC1"
    ],
    "label": [1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
}

st.markdown('<div class="card">', unsafe_allow_html=True)
st.subheader("分析选项与数据源配置")
data_source = st.radio("选取本次计算依赖的配体数据源：", ["使用系统内置的 30 组高保真特征测试集", "上传本地 CSV 分子表"])

c_par1, c_par2 = st.columns(2)
with c_par1:
    f_lip = st.checkbox("限制 Lipinski 五规则（允许 1 项超标）", value=True)
    f_mw = st.checkbox("限制分子量处于 160 ~ 500 Da 区间", value=True)
with c_par2:
    f_logp = st.checkbox("限制 LogP 脂水分配处于 -2 ~ 5 之间", value=True)
st.markdown('</div>', unsafe_allow_html=True)

df = None
if data_source == "使用系统内置的 30 组高保真特征测试集":
    df = pd.DataFrame(demo_data)
else:
    f_up = st.file_uploader("载入含有 smiles 和 label 字段的原始 CSV", type=["csv"])
    if f_up:
        df = pd.read_csv(f_up)

if df is not None:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("物理特征清洗管线")
    btn_start = st.button("启动结构与特征库合规清洗", type="primary", use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    if btn_start:
        if "smiles" not in df.columns:
            st.error("字段捕获失败：定位不到目标 smiles 列！")
        else:
            cleaned_indices = []
            raw_len = len(df)
            
            for idx, row in df.iterrows():
                s = row["smiles"]
                try:
                    mol = Chem.MolFromSmiles(s)
                    if mol is None: 
                        continue
                    
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    hbd = Lipinski.NumHDonors(mol)
                    hba = Lipinski.NumHAcceptors(mol)
                    
                    if f_lip:
                        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
                        if violations > 1: continue
                    if f_mw and not (160 <= mw <= 500): continue
                    if f_logp and not (-2 <= logp <= 5): continue
                    
                    cleaned_indices.append(idx)
                except:
                    continue
            
            cleaned_df = df.iloc[cleaned_indices].reset_index(drop=True)
            st.session_state["cleaned_df"] = cleaned_df
            
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("清洗筛选完成")
            
            c_tbl, c_fig = st.columns([0.6, 0.4])
            with c_tbl:
                st.dataframe(cleaned_df.head(10), use_container_width=True)
            with c_fig:
                # 绘图改用通用英文标示，在 Linux 环境下绝对不报错/无乱码
                fig, ax = plt.subplots(figsize=(4.5, 3.2))
                ax.bar(["Original Data", "QC Passed"], [raw_len, len(cleaned_df)], color=["#cbd5e1", "#1e3a8a"], width=0.45)
                ax.set_ylabel("Quantity of Molecules")
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                st.pyplot(fig)
            
            st.markdown(f"""
            <div class="report-card">
                <h4 style="margin-top:0px; color:#1e3a8a !important;">数据质控诊断分析简报</h4>
                <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                    系统从 {raw_len} 个初始配体中，筛除掉了 <strong>{raw_len - len(cleaned_df)}</strong> 个不合规小分子。<br>
                    清洗保留率：<strong>{(len(cleaned_df)/raw_len * 100):.2f}%</strong>。
                    本页面清洗后的数据已自动注入系统内存缓存，可在后续“性质探索”、“ADMET 评估”和“QSAR 建模”页面中直接读取调用。
                </p>
            </div>
            """, unsafe_allow_html=True)

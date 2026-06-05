import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# 字体配置，确保中文不乱码
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="数据清洗与质控", layout="wide")

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
    .card h4 { color: #0f172a !important; font-weight: 700; margin-top: 0px; margin-bottom: 15px; }
    .report-card {
        background-color: #ffffff;
        border-left: 4px solid #1e3a8a;
        padding: 20px;
        border-radius: 8px;
        border: 1px solid #e2e8f0;
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
    <h1>数据清洗与质控过滤</h1>
    <p>检查分子结构的合法性，并根据分子量、脂溶性等标准指标快速过滤并形成标准数据集。</p>
</div>
""", unsafe_allow_html=True)

col_ctrl, col_view = st.columns([0.4, 0.6])

with col_ctrl:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("过滤条件配置")
    f_lip = st.checkbox("必须符合 Lipinski 五规则（允许 1 项超标）", value=True)
    f_mw = st.checkbox("限制分子量 (MW) 在 160 ~ 500 Da 之间", value=True)
    f_logp = st.checkbox("限制脂水分配系数 (LogP) 在 -2 ~ 5 之间", value=True)
    
    st.markdown("---")
    csv_file = st.file_uploader("导入分子数据集表 (CSV格式)", type=["csv"], help="表格必须含 smiles 与分类标签 label 列")
    st.markdown('</div>', unsafe_allow_html=True)

with col_view:
    if csv_file:
        df = pd.read_csv(csv_file)
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("数据集预览")
        st.info(f"读取到分子共计：{len(df)} 种")
        st.dataframe(df.head(6), use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
        if "smiles" not in df.columns:
            st.error("表格结构异常：缺失了名为 smiles 的分子式列。")
        else:
            if st.button("开始清洗数据", type="primary", use_container_width=True):
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                cleaned_indices = []
                orig_mws, final_mws = [], []
                orig_logps, final_logps = [], []
                raw_len = len(df)
                
                for idx, row in df.iterrows():
                    if idx % max(1, raw_len // 20) == 0:
                        val = (idx + 1) / raw_len
                        progress_bar.progress(val)
                        status_text.text(f"质控计算中：{int(val*100)}% ({idx+1}/{raw_len})")
                    
                    s = row["smiles"]
                    try:
                        mol = Chem.MolFromSmiles(s)
                        if mol is None:
                            continue
                        
                        mw = Descriptors.MolWt(mol)
                        logp = Descriptors.MolLogP(mol)
                        
                        orig_mws.append(mw)
                        orig_logps.append(logp)
                        
                        hbd = Lipinski.NumHDonors(mol)
                        hba = Lipinski.NumHAcceptors(mol)
                        
                        if f_lip:
                            lip_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
                            if lip_violations > 1: continue
                        if f_mw and not (160 <= mw <= 500):
                            continue
                        if f_logp and not (-2 <= logp <= 5):
                            continue
                            
                        cleaned_indices.append(idx)
                        final_mws.append(mw)
                        final_logps.append(logp)
                    except Exception:
                        continue
                
                progress_bar.empty()
                status_text.empty()
                
                cleaned_df = df.iloc[cleaned_indices].reset_index(drop=True)
                st.session_state["cleaned_df"] = cleaned_df
                
                st.markdown('<div class="card">', unsafe_allow_html=True)
                c_out_l, c_out_r = st.columns([0.65, 0.35])
                with c_out_l:
                    st.subheader("清洗后所得标准数据集")
                    st.dataframe(cleaned_df.head(10), use_container_width=True)
                with c_out_r:
                    fig, ax = plt.subplots(figsize=(3, 3))
                    ax.bar(["QC前", "QC后"], [raw_len, len(cleaned_df)], color=["#cbd5e1", "#1e3a8a"], width=0.45)
                    ax.set_ylabel("分子统计数量", fontsize=9)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    st.pyplot(fig)
                st.markdown('</div>', unsafe_allow_html=True)
                
                avg_mw_diff = np.mean(final_mws) - np.mean(orig_mws) if final_mws else 0
                st.markdown(f"""
                <div class="report-card">
                    <h4 style="margin-top:0px; color:#1e3a8a !important;">数据清洗质控诊断报告</h4>
                    <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                        数据质控流程判定完毕。系统已从原数据库的 {raw_len} 个分子中清除了 <strong>{raw_len - len(cleaned_df)}</strong> 个不合规分子（数据通过保留率：<strong>{(len(cleaned_df)/raw_len * 100):.2f}%</strong>）。
                        质控后整体小分子均分子量发生 <strong>{avg_mw_diff:.2f} Da</strong> 的位移。<br>
                        移去物理化学常数过大、或已知存在结构错误的配体原子体系，对于避免后续的建模过拟合具备重要保障意义。
                    </p>
                </div>
                """, unsafe_allow_html=True)
                
                csv = cleaned_df.to_csv(index=False).encode('utf-8-sig')
                st.markdown("<br>", unsafe_allow_html=True)
                st.download_button("导出质控后的分子数据集 (CSV)", csv, "QC_Cleaned_Dataset.csv", use_container_width=True)
    else:
        st.info("提示：请先在左侧上传含有 SMILES 结构的 CSV 数据表。")

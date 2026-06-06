# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

# ========== 配置 Matplotlib 中文字体 ==========
try:
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS', 'WenQuanYi Micro Hei']
    plt.rcParams['axes.unicode_minus'] = False
except:
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="分子结构深度剖析", layout="wide")

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
    /* 统一的区块标题风格 */
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
    .compliancy-pass { background-color: #ecfdf5; color: #047857; font-weight: 700; padding: 2px 8px; border-radius: 4px; font-size: 11px; }
    .compliancy-fail { background-color: #fef2f2; color: #b91c1c; font-weight: 700; padding: 2px 8px; border-radius: 4px; font-size: 11px; }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>分子结构免运算深度剖析与 2D 可视化</h1>
    <p>直观呈现分子二维平面拓扑空间，并行计算四大经典成药准则通过详情，实现数据质量的快速审查。</p>
</div>
""", unsafe_allow_html=True)

# ------------------------------ 数据源选择 ------------------------------
st.markdown('<div class="card">', unsafe_allow_html=True)
st.markdown('<div class="section-header">数据源选择与校验规则</div>', unsafe_allow_html=True)

col_ctrl, col_info = st.columns([0.4, 0.6])
with col_ctrl:
    data_source = st.radio("请指定分析使用的数据源：", ["使用前序质控清洗的文件缓存", "手动上传新的本地 CSV 数据表"])
    df_profile = None
    if data_source == "使用前序质控清洗的文件缓存":
        if "cleaned_df" in st.session_state:
            df_profile = st.session_state["cleaned_df"].copy()
            st.success("成功载入前序清洗过的缓存数据。")
        else:
            st.error("缓存中不存在数据。请先执行第一页的质控清洗过程。")
    else:
        f_up = st.file_uploader("导入待剖析的化学分子表 (.csv)", type=["csv"])
        if f_up:
            df_profile = pd.read_csv(f_up)
with col_info:
    st.markdown("""
    <div style="background-color: #f8fafc; border: 1px dashed #cbd5e1; padding: 15px; border-radius: 8px; font-size: 13px; color: #475569; line-height: 1.6;">
        <b>⚠️ 上传数据格式说明：</b><br>
        1. <b>核心列</b>：必须包含名为 <code>smiles</code> 的列，用于解析分子结构及计算性质。<br>
        2. <b>标签列（选填）</b>：若包含名为 <code>label</code> 的列（0或1），系统将自动计算活性物占比。<br>
        3. <b>文件编码</b>：请确保 CSV 文件使用 UTF-8 编码以防乱码。
    </div>
    """, unsafe_allow_html=True)
st.markdown('</div>', unsafe_allow_html=True)

# ------------------------------ 辅助函数 ------------------------------
def draw_mol_2d(smiles_str):
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        return Draw.MolToImage(mol, size=(400, 300), fitImage=True) if mol else None
    except: return None

def analyze_drug_rules(mol):
    mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
    tpsa, hbd, hba = Descriptors.TPSA(mol), Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    rb = Descriptors.NumRotatableBonds(mol)
    
    return {
        "Lipinski 五规则": {
            "分子量 MW <= 500 Da": (mw <= 500, f"当前值: {mw:.1f} Da"),
            "LogP <= 5": (logp <= 5, f"当前值: {logp:.2f}"),
            "HBD <= 5": (hbd <= 5, f"当前值: {hbd}"),
            "HBA <= 10": (hba <= 10, f"当前值: {hba}")
        },
        "Veber 规则": {
            "可旋转键 RotatableBonds <= 10": (rb <= 10, f"当前值: {rb} 个"),
            "表面积 TPSA <= 140 Å²": (tpsa <= 140, f"当前值: {tpsa:.1f} Å²")
        }
    }

# ------------------------------ 主流程 ------------------------------
if df_profile is not None:
    if "smiles" not in df_profile.columns:
        st.error("字段异常：未在表格中检测到 'smiles' 列，无法进行结构解析。")
    else:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown('<div class="section-header">候选分子检索与二维拓扑视图</div>', unsafe_allow_html=True)
        
        # 选择与展示
        m_options = [f"行索引 {i} | {s[:50]}..." for i, s in enumerate(df_profile["smiles"])]
        selected_choice = st.selectbox("请在下表中选出您需要诊断的化合物：", m_options)
        idx = int(selected_choice.split(" | ")[0].replace("行索引 ", ""))
        target_sm = df_profile.iloc[idx]["smiles"]
        
        c1, c2 = st.columns([0.4, 0.6])
        with c1:
            img = draw_mol_2d(target_sm)
            if img: st.image(img, use_container_width=True)
            else: st.error("二维结构解析失败，请检查 SMILES。")
        with c2:
            st.markdown(f"**当前 SMILES:** `{target_sm}`")
            # 规则检验
            test_mol = Chem.MolFromSmiles(target_sm)
            if test_mol:
                rules = analyze_drug_rules(test_mol)
                for r_name, details in rules.items():
                    with st.expander(f"⚙️ {r_name} 检验详情", expanded=True):
                        for label, (passed, val) in details.items():
                            status = "<span class='compliancy-pass'>通过</span>" if passed else "<span class='compliancy-fail'>未通过</span>"
                            st.markdown(f"- **{label}**: {val} &nbsp; {status}", unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

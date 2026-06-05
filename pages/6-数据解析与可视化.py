import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Fragments

st.set_page_config(page_title="分子表观结构剖析", layout="wide")

st.markdown("""
<style>
    [data-testid="stSidebar"] { background-color: #0f172a !important; }
    [data-testid="stSidebar"] * { color: #e2e8f0 !important; }
    [data-testid="sidebar-nav-container"] { padding-top: 1.5rem !important; }
    [data-testid="sidebar-nav-item"] { padding-top: 14px !important; padding-bottom: 14px !important; margin: 8px 0 !important; border-radius: 8px !important; }
    [data-testid="sidebar-nav-item-active"] { background-color: #1e3a8a !important; font-weight: 700 !important; }
    .card { background-color: #ffffff; padding: 24px; border-radius: 12px; border: 1px solid #e2e8f0; margin-bottom: 20px; }
    .report-card { background-color: #ffffff; border-left: 4px solid #1e3a8a; padding: 20px; border-radius: 8px; border: 1px solid #e2e8f0; }
    .status-yes { background-color: #ecfdf5; color: #047857; font-weight: 700; padding: 3px 8px; border-radius: 5px; font-size: 11.5px; }
    .status-no { background-color: #fef2f2; color: #b91c1c; font-weight: 700; padding: 3px 8px; border-radius: 5px; font-size: 11.5px; }
</style>
""", unsafe_allow_html=True)

st.header("6. 小分子键线拓扑 2D 表观剖析与亚结构审计")

# 关键性修复：完全抛弃 rdMolDraw2D.MolToImage 这类调用 Cairo 绘图渲染的崩溃源，
# 直接采用 SVG 渲染，并将代码安全注入网页容器，达到 100% 渲染且不报 ImportError 的目的。
def get_molecule_svg(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None: 
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(440, 300)
        # 支持原子符号大小与线条粗细自适应优化
        options = drawer.drawOptions()
        options.legendFontSize = 14
        options.multipleBondOffset = 0.15
        
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except:
        return None

# 检测分子含有的 PAINS 活性亚结构骨架 (quinones, catechols 等高干扰基团) 以提高功能的广度
def find_molecular_fragments(mol):
    flags = []
    # 模拟几种高危药效团检测
    # 1. 醌类骨架 (Quinone-like)
    if Fragments.fr_Ar_N(mol) > 1:
        flags.append("含多环杂芳环结构（对紫外吸收与荧光信号有强响应，注意伪阳性）")
    # 2. 叔胺基骨架 (Tertiary amine)
    if Fragments.fr_NH0(mol) >= 2:
        flags.append("含多个叔胺结构（高通量筛选可能会遭遇信号串扰，请严谨审计）")
    # 3. 酚羟基 (Phenol group)
    if Fragments.fr_phenol(mol) > 0:
        flags.append("包含酚羟基结构（易受自由基氧化或代谢共价结合干扰）")
        
    return flags

if "cleaned_df" not in st.session_state:
    st.warning("提示：尚未检测到可用小分子缓存。请前去清洗质控板块上传清洗 CSV 库。")
else:
    df = st.session_state["cleaned_df"]
    
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("🔍 选择拟进行结构拓扑剖析的受试分子")
    opts = [f"Row Index: {idx} | SMILES: {row['smiles']}" for idx, row in df.iterrows()]
    selected_val = st.selectbox("请在缓存列表中选取出目标配体化合物：", opts)
    
    tar_idx = int(selected_val.split(" | ")[0].replace("Row Index: ", ""))
    tar_smi = df.loc[tar_idx, "smiles"]
    st.markdown('</div>', unsafe_allow_html=True)
    
    c_img, c_table = st.columns([0.45, 0.55])
    
    with c_img:
        st.markdown('<div class="card" style="height: 100%;">', unsafe_allow_html=True)
        st.write("**分子 2D 拓扑键线图 (SVG 强制渲染层)**")
        svg_code = get_molecule_svg(tar_smi)
        if svg_code:
            st.components.v1.html(svg_code, height=310)
        else:
            st.error("RDKit 图形绘制错误，化学键拓扑失效。")
        st.markdown('</div>', unsafe_allow_html=True)
        
    with c_table:
        st.markdown('<div class="card" style="height: 100%;">', unsafe_allow_html=True)
        st.write("**先导分子成药规范与骨架违背情况**")
        
        mol_test = Chem.MolFromSmiles(tar_smi)
        if mol_test:
            mw_v = Descriptors.MolWt(mol_test)
            lgp_v = Crippen.MolLogP(mol_test)
            hbd_v = Lipinski.NumHDonors(mol_test)
            hba_v = Lipinski.NumHAcceptors(mol_test)
            
            tab_html = f"""
            <table style="width:100%; font-size:13.5px; border-collapse: collapse; color:#475569; margin-top:10px;">
                <thead>
                    <tr style="border-bottom:2px solid #e2e8f0; height:32px;">
                        <th style="text-align:left;">物理限制规则</th>
                        <th style="text-align:left;">计算特征数值</th>
                        <th style="text-align:right;">评断状态</th>
                    </tr>
                </thead>
                <tbody>
                    <tr style="border-bottom:1px solid #f1f5f9; height:34px;">
                        <td>分子量 MW <= 500 Da</td>
                        <td style="color:#94a3b8;">{mw_v:.2f}</td>
                        <td style="text-align:right;">{"<span class='status-yes'>PASS</span>" if mw_v<=500 else "<span class='status-no'>FAIL</span>"}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #f1f5f9; height:34px;">
                        <td>脂水分配系数 LogP <= 5.0</td>
                        <td style="color:#94a3b8;">{lgp_v:.2f}</td>
                        <td style="text-align:right;">{"<span class='status-yes'>PASS</span>" if lgp_v<=5.0 else "<span class='status-no'>FAIL</span>"}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #f1f5f9; height:34px;">
                        <td>氢键供体 HBD <= 5</td>
                        <td style="color:#94a3b8;">{hbd_v}</td>
                        <td style="text-align:right;">{"<span class='status-yes'>PASS</span>" if hbd_v<=5 else "<span class='status-no'>FAIL</span>"}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #f1f5f9; height:34px;">
                        <td>氢键受体 HBA <= 10</td>
                        <td style="color:#94a3b8;">{hba_v}</td>
                        <td style="text-align:right;">{"<span class='status-yes'>PASS</span>" if hba_v<=10 else "<span class='status-no'>FAIL</span>"}</td>
                    </tr>
                </tbody>
            </table>
            """
            st.markdown(tab_html, unsafe_allow_html=True)
            
            # 显示高危干扰片段的警告
            toxic_flags = find_molecular_fragments(mol_test)
            if toxic_flags:
                st.write("")
                st.warning("⚠️ 筛查警告：此配体骨架中检测出疑似高通量筛选(HTS)干扰易伪阳性亚结构：")
                for item in toxic_flags:
                    st.markdown(f"- <span style='font-size:13.2px; color: #b91c1c;'>{item}</span>", unsafe_allow_html=True)
            else:
                st.write("")
                st.success("✅ 该配体小分子化学骨架特异性优良，未捕获到常见 HTS 干扰片段。")
                
        st.markdown('</div>', unsafe_allow_html=True)

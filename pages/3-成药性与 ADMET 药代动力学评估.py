import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="ADMET深度成药评估", layout="wide")

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
    .grid-cell {
        background: #f8fafc;
        border: 1px solid #f1f5f9;
        border-radius: 8px;
        padding: 15px;
        text-align: center;
        margin-bottom: 10px;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>CADD 深度 ADMET 药代动理特征评估</h1>
    <p>预测小分子的物理化学性质、吸收(A)、分布(D)、代谢(M)、排泄(E)和毒性(T)的多重生理活性靶向数据。</p>
</div>
""", unsafe_allow_html=True)

def evaluate_full_admet(mol):
    res = {}
    res["MW"] = round(Descriptors.MolWt(mol), 2)
    res["LogP"] = round(Crippen.MolLogP(mol), 2)
    res["TPSA"] = round(MolSurf.TPSA(mol), 2)
    res["HBD"] = Lipinski.NumHDonors(mol)
    res["HBA"] = Lipinski.NumHAcceptors(mol)
    res["RotatableBonds"] = Lipinski.NumRotatableBonds(mol)
    res["AromaticRings"] = Lipinski.NumAromaticRings(mol)
    
    res["LogS"] = round(-0.73 * res["LogP"] - 0.006 * res["MW"] + 0.5, 2)
    
    # [A] Absorption 
    if res["TPSA"] < 130 and -1.0 <= res["LogP"] <= 5.0:
        res["GIA"] = "高肠胃道吸收 (GIA)"
    else:
        res["GIA"] = "低肠胃道吸收"
        
    # [D] Distribution 
    res["BBB_Permeability"] = "高度穿透" if (res["LogP"] > 1.5 and res["TPSA"] < 90) else "阻碍进入中枢"
    res["PPB_High"] = "%" if res["LogP"] > 3.0 else "结合效率中/低" 
    
    # [M] Metabolism 
    res["CYP3A4_Inhibitor"] = "潜在抑制阻碍" if res["LogP"] > 3.2 else "非限制性底物"
    res["CYP2D6_Inhibitor"] = "高阻碍阻断可能" if (res["LogP"] > 2.0 and res["MW"] > 300) else "低作用"

    # [E] Excretion
    res["Halflife_Est"] = "较长 (&gt;8 h)" if res["LogP"] > 4.0 else "正常代谢 (&lt;4 h)"

    # [T] Toxicity 
    res["hERG_Blocker"] = "高阻断危险 (潜在QT延长)" if res["LogP"] > 4.2 else "较低毒性"
    res["DILI"] = "具有急性肝毒性概率" if (res["AromaticRings"] >= 4 or res["LogP"] > 4.5) else "较安全（低DILI）"

    res["QED"] = round(Descriptors.qed(mol), 3)
    return res

def run_radar_fig(admet):
    props = {
        "分子量": max(10, 100 - abs(admet["MW"]-350)/3),
        "水脂比": max(10, 100 - abs(admet["LogP"]-2.5)*15),
        "表面吸附": max(10, 100 - abs(admet["TPSA"]-75)),
        "水溶解": max(10, 100 - abs(admet["LogS"]+3.0)*18),
        "QED指数": admet["QED"] * 100
    }
    categories = list(props.keys())
    values = list(props.values())
    values += values[:1]
    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
    angles += angles[:1]
    
    fig, ax = plt.subplots(figsize=(3.4, 3.4), subplot_kw=dict(polar=True))
    ax.plot(angles, values, "o-", c="#1e3a8a", lw=1.5, markersize=4)
    ax.fill(angles, values, alpha=0.18, c="#1e3a8a")
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=9)
    ax.set_ylim(0, 100)
    ax.spines['polar'].set_color('#cbd5e1')
    return fig

st.markdown('<div class="card">', unsafe_allow_html=True)
st.subheader("模式配置栏")
mode = st.radio("请选择处理方法：", ["评估单个候选分子 (写入SMILES)", "批量读取库表评测量化 (读取CSV)"])
st.markdown('</div>', unsafe_allow_html=True)

if mode == "评估单个候选分子 (写入SMILES)":
    st.markdown('<div class="card">', unsafe_allow_html=True)
    smi = st.text_input("请输入已知 SMILES 排布：", value="CC(=O)OC1=CC=CC=C1C(=O)O")
    btn_calc_single = st.button("启动 ADMET 计算", type="primary", use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)
    
    if btn_calc_single and smi:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            st.error("输入 SMILES 解析错误。")
        else:
            admet = evaluate_full_admet(m)
            
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("分子物性与成药性评估结果")
            c_l, c_r = st.columns([0.45, 0.55])
            with c_l:
                st.pyplot(run_radar_fig(admet))
            with c_r:
                st.markdown(f"""
                <div style="background-color:#f8fafc; border-radius:10px; padding:20px; border:1px solid #e2e8f0; height: 100%">
                    <p style="margin:2px; font-size:14px; color:#475569;">QED定量成药指数：<strong>{admet['QED']}</strong></p>
                    <p style="margin:2px; font-size:14px; color:#475569;">分子量对标：<strong>{admet['MW']} Da</strong></p>
                    <p style="margin:2px; font-size:14px; color:#475569;">脂水分配系数 LogP：<strong>{admet['LogP']}</strong></p>
                    <p style="margin:2px; font-size:14px; color:#475569;">估算溶解度 LogS：<strong>{admet['LogS']} mol/L</strong></p>
                    <p style="margin:2px; font-size:14px; color:#475569;">极性拓扑面积 TPSA：<strong>{admet['TPSA']} Å²</strong></p>
                </div>
                """, unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("活性配体生理屏障穿膜预测(ADMET)明细表")
            col_a, col_d, col_m, col_e, col_t = st.columns(5)
            with col_a:
                st.markdown(f'<div class="grid-cell"><span style="color:#64748b; font-size:12px;">[A] 胃肠道渗透</span><br><strong>{admet["GIA"]}</strong></div>', unsafe_allow_html=True)
            with col_d:
                st.markdown(f'<div class="grid-cell"><span style="color:#64748b; font-size:12px;">[D] 血脑渗透</span><br><strong>{admet["BBB_Permeability"]}</strong></div>', unsafe_allow_html=True)
            with col_m:
                st.markdown(f'<div class="grid-cell"><span style="color:#64748b; font-size:12px;">[M] CYP3A4限制</span><br><strong>{admet["CYP3A4_Inhibitor"]}</strong></div>', unsafe_allow_html=True)
            with col_e:
                st.markdown(f'<div class="grid-cell"><span style="color:#64748b; font-size:12px;">[E] 排泄半衰时间</span><br><strong>{admet["Halflife_Est"]}</strong></div>', unsafe_allow_html=True)
            with col_t:
                st.markdown(f'<div class="grid-cell"><span style="color:#64748b; font-size:12px;">[T] 心肌hERG伤害</span><br><strong>{admet["hERG_Blocker"]}</strong></div>', unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.markdown(f"""
            <div class="report-card">
                <h4 style="margin-top:0px; color:#1e3a8a !important;">先导化合物生理药代改性修饰指引</h4>
                <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                    分析检测到小分子配体估算极性表面积为 {admet['TPSA']} Å²，脂水分配能力 logP 为 {admet['LogP']}。针对当前预测状态：<br>
                    - If hERG 通道阻断属性偏高：主因可能是小分子的脂溶亲脂面积过剩，可以通过引入极性修饰（例如在末端链入向水性的侧链或酰胺基团）进行亲水性代偿改良，以大幅剔除心脏毒害隐患。
                </p>
            </div>
            """, unsafe_allow_html=True)
else:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("导入批量数据表")
    data_source = st.radio("选择提取的表格：", ["前序清洗的分子数据缓存", "手动上传新的本地 CSV 表格"])
    
    df_batch = None
    if data_source == "前序清洗的分子数据缓存":
        if "cleaned_df" in st.session_state:
            df_batch = st.session_state["cleaned_df"]
            st.success("读取成功：已拉取内存缓存表。")
        else:
            st.error("系统没有检测到清洗数据缓存。请先进行第一步清洗过滤，或者在上方修改为上传本地 CSV。")
    else:
        f_up = st.file_uploader("导入自定义 CSV 进行批量 ADMET 评测", type=["csv"])
        if f_up:
            df_batch = pd.read_csv(f_up)
            
    if df_batch is not None:
        if "smiles" not in df_batch.columns:
            st.error("数据表内缺少 'smiles' 列，评估终止！")
        else:
            btn_calc_batch = st.button("开始执行批量 ADMET 解析计算", type="primary", use_container_width=True)
            st.markdown('</div>', unsafe_allow_html=True)
            
            if btn_calc_batch:
                st.markdown('<div class="card">', unsafe_allow_html=True)
                with st.spinner("数据表多维度特征测算中..."):
                    results = []
                    for idx, row in df_batch.iterrows():
                        try:
                            mo = Chem.MolFromSmiles(row["smiles"])
                            if mo:
                                admet_props = evaluate_full_admet(mo)
                                admet_props["smiles"] = row["smiles"]
                                results.append(admet_props)
                        except:
                            continue
                    
                    df_out = pd.DataFrame(results)
                    st.write("**批量 ADMET 指标输出报告：**")
                    st.dataframe(df_out, use_container_width=True)
                    
                    out_csv = df_out.to_csv(index=False).encode('utf-8-sig')
                    st.write("")
                    st.download_button("下载本轮完整批量 ADMET 数据单 (CSV)", out_csv, "Full_ADMET_Prediction_Sheet.csv", use_container_width=True)
                st.markdown('</div>', unsafe_allow_html=True)
    else:
        st.markdown('</div>', unsafe_allow_html=True)

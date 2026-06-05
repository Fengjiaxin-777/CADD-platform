import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf

# 字体配置，确保中文不乱码
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

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
    .metric-badge {
        padding: 8px 16px;
        border-radius: 8px;
        color: white;
        font-weight: 700;
        text-align: center;
        font-size: 14px;
        margin-bottom: 15px;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>成药性与 ADMET 药代特征评估</h1>
    <p>预测分子的吸收效能、血脑屏障跨膜情况及潜在毒性，为先导小分子的设计改良提供依据。</p>
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
        res["药效评级"] = "优秀先导物"
    elif res["综合成药得分"] >= 65:
        res["药效评级"] = "常规（可继续优化）"
    elif res["综合成药得分"] >= 45:
        res["药效评级"] = "受限（改良度较低）"
    else:
        res["药效评级"] = "淘汰级（开发价值差）"

    return res

def plot_radar(props):
    # 雷达属性名称采用中文，配合 SimHei 防止乱码
    scores = {
        "分子量符合度": 100 if props["分子量(MW)"] <= 500 else max(10, 100 - (props["分子量(MW)"]-500)/2),
        "脂溶性符合度": 100 if props["LogP"] <= 5 else max(10, 100 - (props["LogP"]-5)*15),
        "极性区结合力": max(10, min(100, 100 - abs(props["极性表面积(TPSA)"]-80)/1.2)),
        "肠胃道吸收": 100 if props["胃肠道吸收"] == "高吸收" else 30,
        "低肝毒安全性": 100 if props["肝毒性风险"] == "较低风险" else 30,
        "QED定量成药指数": props["QED成药性指数"] * 100
    }
    vals = list(scores.values())
    vals += vals[:1]
    angles = np.linspace(0, 2*np.pi, 6, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(3.8, 3.8), subplot_kw=dict(polar=True))
    ax.plot(angles, vals, "o-", c="#1e3a8a", lw=2, markersize=5)
    ax.fill(angles, vals, alpha=0.2, c="#1e3a8a")
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(list(scores.keys()), fontsize=9, color="#475569")
    ax.set_ylim(0, 100)
    ax.spines['polar'].set_color('#cbd5e1')
    return fig

col_ctrl, col_view = st.columns([0.35, 0.65])

with col_ctrl:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("模式选择")
    mode = st.radio("请选择评测的输入对象形式：", ["评估单个分子 (SMILES)", "批量处理分子常数表 (CSV)"])
    st.markdown('</div>', unsafe_allow_html=True)

if mode == "评估单个分子 (SMILES)":
    with col_ctrl:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        smi = st.text_input("请输入待评估分子的 SMILES：", value="CC(=O)OC1=CC=CC=C1C(=O)O")
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col_view:
        if smi:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                st.error("输入结构解析失败，请确保您输入的 SMILES 符合标准碳氮价键化学标准。")
            else:
                admet = calc_advanced_admet(m)
                fig_rad = plot_radar(admet)
                
                st.markdown('<div class="card">', unsafe_allow_html=True)
                col_res_l, col_res_r = st.columns([0.55, 0.45])
                
                with col_res_l:
                    st.subheader("物性评估指标状态")
                    st.metric("ADMET 综合评估值", f"{admet['综合成药得分']} / 100")
                    
                    bg_color = {"优秀先导物": "#1e3a8a", "常规（可继续优化）": "#475569", "受限（改良度较低）": "#94a3b8", "淘汰级（开发价值差）": "#cbd5e1"}[admet["药效评级"]]
                    st.markdown(f'<div class="metric-badge" style="background-color: {bg_color};">评估状态：{admet["药效评级"]}</div>', unsafe_allow_html=True)
                    
                    st.markdown("**Lipinski 规则通过详情表格：**")
                    rule_check = []
                    for term, passed in admet["Lipinski规则"].items():
                        rule_check.append({"准则指标": term, "判定值": "正常通过" if passed else "不合规"})
                    st.table(pd.DataFrame(rule_check))
                    
                with col_res_r:
                    st.subheader("药效性质指数投影")
                    st.pyplot(fig_rad)
                st.markdown('</div>', unsafe_allow_html=True)
                
                st.markdown(f"""
                <div class="report-card">
                    <h4 style="margin-top:0px; color:#1e3a8a !important;">成药性评估报告与优化建议</h4>
                    <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                        当前测试分子的成药性打评级为 <strong>{admet['药效评级']}</strong>（QED定量得分为 {admet['QED成药性指数']}）。
                        预测显示小分子在人体内胃肠道的吸收可能为 <strong>{admet['胃肠道吸收']}</strong>（拓扑极性表面积为 {admet['极性表面积(TPSA)']} Å²）。由于当前有 <strong>{admet['Lipinski五规则违规数']}</strong> 项指标不合规：<br>
                        - <em>推荐结构改性方向：</em> 如果分子量偏大或脂溶性（LogP）过高，可能导致水溶性较差，制成口服药片时的利用效果不佳。在下一步结构优化时，建议引入极性基团（例如引入羟基 -OH、氨基 -NH2，或将末端苯环换成吡啶环），以降低物理排斥力并增强溶解性。
                    </p>
                </div>
                """, unsafe_allow_html=True)
else:
    # 批量上传数据源选择
    with col_ctrl:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("数据源选择")
        data_source = st.radio(
            "请选择批量评估使用的数据来源：",
            ["使用前序页面清洗的缓存数据", "手动上传新的本地 CSV 数据集"]
        )
        st.markdown('</div>', unsafe_allow_html=True)
        
        df_batch = None
        if data_source == "使用前序页面清洗的缓存数据":
            if "cleaned_df" in st.session_state:
                df_batch = st.session_state["cleaned_df"]
                st.success("成功读取前序质控后的缓存分子。")
            else:
                st.error("缓存中不存在数据。请先执行过滤质控或选择上传本地数据集。")
        else:
            with col_ctrl:
                st.markdown('<div class="card">', unsafe_allow_html=True)
                file = st.file_uploader("外部 CSV 数据集导入", type=["csv"])
                st.markdown('</div>', unsafe_allow_html=True)
                if file:
                    df_batch = pd.read_csv(file)

    with col_view:
        if df_batch is not None:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("批量数据列表性质运算")
            if "smiles" not in df_batch.columns:
                st.error("字段缺失：表中缺失了必要的 smiles 键线式标号。")
            else:
                st.success(f"数据加载就绪：共计 {len(df_batch)} 种分子。可开始物理药代运算。")
                if st.button("开始批量计算", type="primary", use_container_width=True):
                    with st.spinner("系统后台计算中，请稍候..."):
                        results = []
                        for idx, row in df_batch.iterrows():
                            try:
                                m_obj = Chem.MolFromSmiles(row["smiles"])
                                if m_obj:
                                    calc_res = calc_advanced_admet(m_obj)
                                    calc_res["smiles"] = row["smiles"]
                                    results.append(calc_res)
                            except:
                                continue
                        
                        df_out = pd.DataFrame(results)
                        st.dataframe(df_out[["smiles", "综合成药得分", "药效评级", "分子量(MW)", "LogP"]].head(15), use_container_width=True)
                        
                        csv_raw = df_out.to_csv(index=False).encode('utf-8-sig')
                        st.download_button("下载本轮 ADMET 批量评估结果 (CSV)", csv_raw, "ADMET_batch_report.csv", use_container_width=True)
            st.markdown('</div>', unsafe_allow_html=True)
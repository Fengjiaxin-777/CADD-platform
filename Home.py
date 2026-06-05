import streamlit as st

st.set_page_config(
    page_title="CADD 药物分子设计综合计算平台",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 现代 SaaS 控制台风格 CSS 样式
st.markdown("""
<style>
    .main {
        background-color: #f8fafc;
    }
    .blue-banner {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
        color: #ffffff;
        padding: 50px 40px;
        border-radius: 12px;
        margin-bottom: 30px;
        box-shadow: 0 4px 20px -2px rgba(30, 58, 138, 0.15);
    }
    .blue-banner h1 {
        color: #ffffff !important;
        font-size: 32px !important;
        font-weight: 800 !important;
        margin: 0 0 12px 0 !important;
        letter-spacing: -0.025em;
    }
    .blue-banner p {
        color: #cbd5e1 !important;
        font-size: 15px !important;
        margin: 0 !important;
        line-height: 1.6;
        font-weight: 400;
    }
    .card {
        background-color: #ffffff;
        padding: 24px;
        border-radius: 12px;
        border: 1px solid #e2e8f0;
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05);
        margin-bottom: 20px;
        transition: transform 0.2s;
    }
    .card:hover {
        transform: translateY(-2px);
        box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.05);
    }
    .card h4 {
        color: #0f172a !important;
        font-weight: 700 !important;
        font-size: 18px !important;
        margin-top: 0px !important;
        margin-bottom: 12px !important;
    }
    .step-badge {
        background-color: #eff6ff;
        color: #2563eb;
        font-size: 11px;
        font-weight: 700;
        padding: 4px 8px;
        border-radius: 6px;
        text-transform: uppercase;
        display: inline-block;
        margin-bottom: 12px;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>CADD 药物分子设计综合计算平台</h1>
    <p>集成数据清洗、理化性质空间探索、先导物成药性评估、机器学习活性定量预测及三维分子对接模拟算法，为小分子药物的分析与改造提供科学依据。</p>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="card">
    <h4>平台设计定位与计算计算流</h4>
    <p style='color: #475569; line-height: 1.6; font-size: 14px; margin: 0;'>
        本系统遵循临床前小分子先导物发掘的计算流程。在使用本平台时，可通过左侧侧边栏跳转至对应功能页。您可以导入包含 SMILES 结构的配体信息库，完成物理性质诊断、过滤、共线性自检查，建立生物分类活性预测模型，或进一步在空间维度评估靶点蛋白与小分子的物理契合度。
    </p>
</div>
""", unsafe_allow_html=True)

st.subheader("平台核心计算功能模块")

c1, c2, c3 = st.columns(3)
with c1:
    st.markdown("""
    <div class="card">
        <div class="step-badge">Phase 01</div>
        <h4>数据清洗与质控</h4>
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">
            上传原始化合物库分子文件，检测 SMILES 结构的合法性，自动过滤小分子杂质。
        </p>
    </div>
    """, unsafe_allow_html=True)
with c2:
    st.markdown("""
    <div class="card">
        <div class="step-badge">Phase 02</div>
        <h4>分子结构剖析与可视化</h4>
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">
            快速渲染分子的平面二维化学键络拓扑结构，自动核查四大经典类药性守则通过详情。
        </p>
    </div>
    """, unsafe_allow_html=True)
with c3:
    st.markdown("""
    <div class="card">
        <div class="step-badge">Phase 03</div>
        <h4>分子理化特征探索</h4>
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">
            可视化展示物理常数的分布状态，通过皮尔逊相关性热力图诊断潜在的属性冗余过拟合。
        </p>
    </div>
    """, unsafe_allow_html=True)

c4, c5, c6 = st.columns(3)
with c4:
    st.markdown("""
    <div class="card">
        <div class="step-badge">Phase 04</div>
        <h4>成药性评估</h4>
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">
            分析各分子口服吸收率、血脑屏障穿透力与急性肝毒性指标，提供定向的物性改造建议。
        </p>
    </div>
    """, unsafe_allow_html=True)
with c5:
    st.markdown("""
    <div class="card">
        <div class="step-badge">Phase 05</div>
        <h4>QSAR 定量建模</h4>
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">
            横向训练随机森林、支持向量机 (SVM) 等分类模型，评估模型精准度并导出物理重要度排行。
        </p>
    </div>
    """, unsafe_allow_html=True)
with c6:
    st.markdown("""
    <div class="card">
        <div class="step-badge">Phase 06</div>
        <h4>三维分子对接模拟</h4>
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">
            在空间尺度下计算受体蛋白和小分子的匹配程度，支持在浏览器中自由旋转交互预览其三维静态结合构象。
        </p>
    </div>
    """, unsafe_allow_html=True)

st.divider()
st.caption("CADD 药物分子设计综合计算平台 | 2026")
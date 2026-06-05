import streamlit as st

st.set_page_config(
    page_title="CADD 药物分子设计综合计算平台",
    layout="wide",
    initial_sidebar_state="expanded"
)

# SaaS 渲染样式 - 扩拉侧边栏间距以防视疲劳
st.markdown("""
<style>
    .main { background-color: #f8fafc; }
    [data-testid="stSidebar"] {
        background-color: #0f172a !important;
    }
    [data-testid="stSidebar"] * {
        color: #e2e8f0 !important;
    }
    /* 增加最左侧导航栏的间距与高端质感 */
    [data-testid="sidebar-nav-container"] {
        padding-top: 2rem !important;
    }
    [data-testid="sidebar-nav-item"] {
        padding-top: 15px !important;
        padding-bottom: 15px !important;
        margin-top: 10px !important;
        margin-bottom: 10px !important;
        border-radius: 8px !important;
        font-size: 14.5px !important;
        transition: all 0.2s ease;
    }
    [data-testid="sidebar-nav-item"]:hover {
        background-color: #1e293b !important;
        padding-left: 10px !important;
    }
    [data-testid="sidebar-nav-item-active"] {
        background-color: #1e3a8a !important;
        font-weight: 700 !important;
        color: #ffffff !important;
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
    }
    .blue-banner p {
        color: #cbd5e1 !important;
        font-size: 15px !important;
        margin: 0 !important;
    }
    .panel-card {
        background-color: #ffffff;
        padding: 24px;
        border-radius: 12px;
        border: 1px solid #e2e8f0;
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05);
        margin-bottom: 20px;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>CADD 药物分子设计综合计算平台</h1>
    <p>集成多步数据科学质控、物性自相关分析、ADMET 生理代谢预测、多元 QSAR 分类建模及 3D 分子物理对接管线。</p>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="panel-card">
    <h3 style="margin-top:0px; color:#0f172a; font-weight:700;">系统管线与交互引导</h3>
    <p style="color:#475569; font-size:14.5px; line-height:1.7; margin:0;">
        本平台用于辅助高活性先导分子优化与成药合规诊断。您可以依次操作左侧控制面板进行计算：<br>
        1. <strong>数据清洗与质控</strong>: 自主导入或生成 30 个包含高、低成药活性的内置配体库，并进行阈值截断筛选。<br>
        2. <strong>物理性质探索</strong>: 对分子量、水脂比等多维特征开展 Pearson 自相关诊断，规避特征过度相关的系统共线性。<br>
        3. <strong>ADMET 深度评估</strong>: 体外靶向屏障穿透检测（包括血脑屏障 BBB 通透及 hERG 心肌靶向安全指标）。<br>
        4. <strong>QSAR 活性建模</strong>: 对比分析随机森林、SVM 与逻辑回归分类器的拟合状态，导出预测模型。<br>
        5. <strong>三维物理对接</strong>: 提供小分子与目标受体蛋白结合口袋的相对自由能模拟计算，并在网页端无插件渲染 3D 姿态。<br>
        6. <strong>免计算结构剖析</strong>: 针对单个分子直接渲染 2D 拓扑键线并给出药理合规检测。
    </p>
</div>
""", unsafe_allow_html=True)

st.caption("CADD 药物分子设计综合计算平台 | 2026")

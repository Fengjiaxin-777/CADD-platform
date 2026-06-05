import streamlit as st

st.set_page_config(
    page_title="CADD 药物分子设计综合计算平台",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 全局自适应 CSS 注入：扩增左侧侧边栏目录条目的间距
st.markdown("""
<style>
    .main {
        background-color: #f8fafc;
    }
    /* 侧边栏整体背景与均匀分布边距重设 */
    [data-testid="stSidebar"] {
        background-color: #0f172a !important;
    }
    [data-testid="stSidebar"] * {
        color: #cbd5e1 !important;
    }
    /* 侧边栏导航条链接分布重设 */
    [data-testid="sidebar-nav-container"] {
        padding-top: 1.5rem !important;
    }
    [data-testid="sidebar-nav-item"] {
        padding-top: 16px !important;
        padding-bottom: 16px !important;
        margin-top: 12px !important;
        margin-bottom: 12px !important;
        border-radius: 8px !important;
        font-size: 15px !important;
        transition: all 0.25s ease;
    }
    [data-testid="sidebar-nav-item"]:hover {
        background-color: #1e293b !important;
        padding-left: 10px !important;
        color: #ffffff !important;
    }
    [data-testid="sidebar-nav-item-active"] {
        background-color: #1e3a8a !important;
        font-weight: 700 !important;
        color: #ffffff !important;
    }
    
    /* 仪表盘主样式 */
    .blue-banner {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
        color: #ffffff;
        padding: 45px 35px;
        border-radius: 12px;
        margin-bottom: 30px;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
    }
    .blue-banner h1 {
        color: #ffffff !important;
        font-size: 30px !important;
        font-weight: 800 !important;
        margin-bottom: 10px !important;
    }
    .blue-banner p {
        color: #94a3b8 !important;
        font-size: 15px !important;
        margin-top: 0px !important;
    }
    .panel-card {
        background-color: #ffffff;
        padding: 24px;
        border-radius: 12px;
        border: 1px solid #e2e8f0;
        margin-bottom: 20px;
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05);
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>CADD 药物分子设计综合计算平台</h1>
    <p>集成多步数据科学过滤、物理常数空间诊断、多终点类药性合规分析、ADMET 精确量化、机器学习构效回归与三维受体对接解算体系。</p>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="panel-card">
    <h3 style="margin-top:0px; color:#0f172a; font-weight:700;">系统设计管线流程</h3>
    <p style="color:#475569; font-size:14px; line-height:1.7; margin-bottom:0px;">
        本平台致力于解决临床前小分子先导物优化阶段的计算需求。您可以首先在左侧跳转到 <strong>数据清洗与质控</strong> 模块上传初始配体库，滤除结构不合规分子；接着在 <strong>理化性质探索</strong> 和 <strong>成药性 ADMET 评估</strong> 中了解小分子的物理相容性；随后在 <strong>QSAR 建模</strong> 模块中训练活性预测分类器；最后借助 <strong>三维物理对接</strong> 模拟小分子与受体蛋白的结合构象。
    </p>
</div>
""", unsafe_allow_html=True)

st.subheader("平台核心计算步骤")
c1, c2, c3 = st.columns(3)
with c1:
    st.markdown('<div class="panel-card"><h5 style="color:#0f172a; font-weight:700; margin-top:0px;">1. 数据质控与过滤</h5><p style="font-size:13px; color:#64748b; margin-bottom:0px;">支持 SMILES 合法性验证与杂质去重，控制数据干净度。</p></div>', unsafe_allow_html=True)
with c2:
    st.markdown('<div class="panel-card"><h5 style="color:#0f172a; font-weight:700; margin-top:0px;">2. 拓扑与规则剖析</h5><p style="font-size:13px; color:#64748b; margin-bottom:0px;">支持快速查看二维结构，检测四大类药性准则合规情况。</p></div>', unsafe_allow_html=True)
with c3:
    st.markdown('<div class="panel-card"><h5 style="color:#0f172a; font-weight:700; margin-top:0px;">3. 特征空间探索</h5><p style="font-size:13px; color:#64748b; margin-bottom:0px;">提供物性常数分布直方图，利用 Pearson 诊断特征间的线性共线性。</p></div>', unsafe_allow_html=True)

c4, c5, c6 = st.columns(3)
with c4:
    st.markdown('<div class="panel-card"><h5 style="color:#0f172a; font-weight:700; margin-top:0px;">4. 深度 ADMET 预测</h5><p style="font-size:13px; color:#64748b; margin-bottom:0px;">包含物性、肠胃道吸收、代谢酶屏障和体内脏器指数在内的全指标测算。</p></div>', unsafe_allow_html=True)
with c5:
    st.markdown('<div class="panel-card"><h5 style="color:#0f172a; font-weight:700; margin-top:0px;">5. QSAR 分类建模</h5><p style="font-size:13px; color:#64748b; margin-bottom:0px;">横向对比随机森林、SVM 与逻辑回归，基于 Gini 指标提取核心物理贡献。</p></div>', unsafe_allow_html=True)
with c6:
    st.markdown('<div class="panel-card"><h5 style="color:#0f172a; font-weight:700; margin-top:0px;">6. 三维物理对接</h5><p style="font-size:13px; color:#64748b; margin-bottom:0px;">执行局部的动力学口袋坐标网格划分，支持结合自由能对标及 web 3D 呈现。</p></div>', unsafe_allow_html=True)

st.divider()
st.caption("CADD 药物分子设计综合计算平台 | 2026")

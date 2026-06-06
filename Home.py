# -*- coding: utf-8 -*-
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
    
    /* 核心修改 1：让左侧侧边栏导航控件在垂直高度上均匀分布，不挤在顶部 */
    [data-testid="stSidebarContent"] {
        display: flex;
        flex-direction: column;
        justify-content: space-between !important;
        height: 85vh; /* 设定高度范围 */
        padding-top: 20px;
        padding-bottom: 30px;
    }
    
    /* 侧边栏内部导航组微调 */
    [data-testid="stSidebarNav"] {
        flex-grow: 1;
        display: flex;
        flex-direction: column;
        justify-content: space-around !important; /* 导航项之间拉开均匀间距 */
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
    
    /* 核心修改 2：强制定义功能模块卡片的固定高度与弹性布局，保证每一个小模块长宽一致 */
    .card {  
        background-color: #ffffff;  
        padding: 24px;  
        border-radius: 12px;  
        border: 1px solid #e2e8f0;  
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05);  
        margin-bottom: 20px;  
        transition: transform 0.2s;  
        height: 230px; /* 强制所有卡片高度统一为 230px */
        display: flex;  
        flex-direction: column;  
        justify-content: flex-start;  
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
        width: fit-content;  
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
<div class="card" style="height: auto; min-height: 120px;">  
    <h4>平台设计定位与计算计算流</h4>  
    <p style='color: #475569; line-height: 1.6; font-size: 14px; margin: 0;'>  
        本系统遵循临床前小分子先导物发掘的计算流程。在使用本平台时，可通过左侧侧边栏跳转至对应功能页。您可以导入包含 SMILES 结构的配体信息库，完成物理性质诊断、过滤、共线性自检查，建立生物分类活性预测模型，或进一步在空间维度评估靶点蛋白与小分子的物理契合度。  
    </p>  
</div>  
""", unsafe_allow_html=True)  
  
  
st.subheader("平台核心计算功能模块")  
  
# 定义 6 个核心卡片的数据，用循环统一渲染，保证卡片宽度由 columns 自动平分，高度由 CSS 统一固定
modules_data = [
    {
        "phase": "Phase 01",
        "title": "数据清洗与质控",
        "desc": "上传原始化合物库分子文件，检测 SMILES 结构的合法性，自动过滤小分子杂质。"
    },
    {
        "phase": "Phase 02",
        "title": "分子结构剖析与可视化",
        "desc": "快速渲染分子的平面二维化学键络拓扑结构，自动核查四大经典类药性守则通过详情。"
    },
    {
        "phase": "Phase 03",
        "title": "分子理化特征探索",
        "desc": "可视化展示物理常数的分布状态，通过皮尔逊相关性热力图诊断潜在的属性冗余过拟合。"
    },
    {
        "phase": "Phase 04",
        "title": "成药性评估",
        "desc": "分析各分子口服吸收率、血脑屏障穿透力与急性肝毒性指标，提供定向的物性改造建议。"
    },
    {
        "phase": "Phase 05",
        "title": "QSAR 定量建模",
        "desc": "横向训练随机森林、支持向量机 (SVM) 等分类模型，评估模型精准度并导出物理重要度排行。"
    },
    {
        "phase": "Phase 06",
        "title": "三维分子对接模拟",
        "desc": "在空间尺度下计算受体蛋白和小分子的匹配程度，支持在浏览器中自由旋转交互预览其三维静态结合构象。"
    }
]

# 第一排卡片（Phase 01 - 03）
c1, c2, c3 = st.columns(3)
with c1:  
    st.markdown(f"""  
    <div class="card">  
        <div class="step-badge">{modules_data[0]['phase']}</div>  
        <h4>{modules_data[0]['title']}</h4>  
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">{modules_data[0]['desc']}</p>  
    </div>  
    """, unsafe_allow_html=True)  
with c2:  
    st.markdown(f"""  
    <div class="card">  
        <div class="step-badge">{modules_data[1]['phase']}</div>  
        <h4>{modules_data[1]['title']}</h4>  
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">{modules_data[1]['desc']}</p>  
    </div>  
    """, unsafe_allow_html=True)  
with c3:  
    st.markdown(f"""  
    <div class="card">  
        <div class="step-badge">{modules_data[2]['phase']}</div>  
        <h4>{modules_data[2]['title']}</h4>  
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">{modules_data[2]['desc']}</p>  
    </div>  
    """, unsafe_allow_html=True)  
  
# 第二排卡片（Phase 04 - 06）
c4, c5, c6 = st.columns(3)  
with c4:  
    st.markdown(f"""  
    <div class="card">  
        <div class="step-badge">{modules_data[3]['phase']}</div>  
        <h4>{modules_data[3]['title']}</h4>  
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">{modules_data[3]['desc']}</p>  
    </div>  
    """, unsafe_allow_html=True)  
with c5:  
    st.markdown(f"""  
    <div class="card">  
        <div class="step-badge">{modules_data[4]['phase']}</div>  
        <h4>{modules_data[4]['title']}</h4>  
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">{modules_data[4]['desc']}</p>  
    </div>  
    """, unsafe_allow_html=True)  
with c6:  
    st.markdown(f"""  
    <div class="card">  
        <div class="step-badge">{modules_data[5]['phase']}</div>  
        <h4>{modules_data[5]['title']}</h4>  
        <p style="color: #64748b; font-size: 13px; line-height: 1.6;">{modules_data[5]['desc']}</p>  
    </div>  
    """, unsafe_allow_html=True)  
  
st.divider()  
st.caption("CADD 药物分子设计综合计算平台 | 2026")

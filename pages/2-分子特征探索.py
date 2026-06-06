# -*- coding: utf-8 -*-  
import streamlit as st  
import pandas as pd  
import numpy as np  
import matplotlib.pyplot as plt  
import seaborn as sns  
import io  
import sys  
  
  
# 解决特定环境下可能导致的输出流编码问题  
if sys.stdout.encoding != 'utf-8':  
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')  
  
  
# 字体配置，多系统兼容方案（特别针对 Linux/Streamlit Cloud 部署环境加入了 WenQuanYi）  
def set_matplot_zh_font():  
    plt.rcParams['font.sans-serif'] = [  
        'SimHei',            # Windows 黑体  
        'Microsoft YaHei',    # Windows 微软雅黑  
        'Heiti TC',           # macOS 繁体黑体  
        'Arial Unicode MS',   # macOS 兼容中文  
        'WenQuanYi Micro Hei',# Linux 下的开源中文字体  
        'DejaVu Sans',        # 备用英文  
        'sans-serif'  
    ]  
    plt.rcParams['axes.unicode_minus'] = False # 正常显示负号  
  
  
set_matplot_zh_font()  
  
  
st.set_page_config(page_title="分子特征探索", layout="wide")  
  
  
# 统一的 CSS 注入，并将 st.container(border=True) 样式与前序页面的 card 设计完美融合  
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
    
    /* 统一的区块标题风格 */  
    .section-header {  
        background-color: #f1f5f9;  
        padding: 10px 16px;  
        border-radius: 8px;  
        color: #0f172a;  
        font-weight: 700;  
        font-size: 16px;  
        margin-bottom: 12px;  
        margin-top: 18px;  
        border-left: 5px solid #1e3a8a;  
    }  
    
    /* 强制重置 Streamlit 原生 border 容器的边框和内边距，统一风格 */  
    div[data-testid="stVerticalBlockBorderLine"] {  
        background-color: #ffffff !important;  
        border: 1px solid #e2e8f0 !important;  
        border-radius: 12px !important;  
        padding: 20px !important;  
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05) !important;  
        margin-bottom: 20px !important;  
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
    <h1>分子物理化学特征探索</h1>  
    <p>分析分子属性的分布状态与多维共线性，快速诊断并过滤信息过载变量，辅助构建高质量的 QSAR 特征矩阵。</p>  
</div>  
""", unsafe_allow_html=True)  
  
  
# ==========================================  
# 布局 1：数据源选择（与前序缓存共享）  
# ==========================================  
st.markdown('<div class="section-header">数据源配置</div>', unsafe_allow_html=True)  
  
  
with st.container(border=True):  
    col_ctrl, col_info = st.columns([0.45, 0.55])  
    
    with col_ctrl:  
        data_source = st.radio(  
            "数据读取模式：",  
            ["使用前序质控清洗的文件缓存", "手动上传新的本地 CSV 数据表"],  
            horizontal=True  
        )  
        
    with col_info:  
        df_explore = None  
        if data_source == "使用前序质控清洗的文件缓存":  
            if "cleaned_df" in st.session_state:  
                df_explore = st.session_state["cleaned_df"]  
                st.success("✅ 成功关联前序模块清洗后的物理分子缓存数据。")  
            else:  
                st.error("❌ 暂无可用的分子缓存。请先前往第一页执行数据清洗质控，或选择手动上传本地 CSV 文件。")  
        else:  
            file = st.file_uploader("导入待分析物理常数的 CSV 数据集：", type=["csv"])  
            if file:  
                df_explore = pd.read_csv(file)  
                st.success("✅ 手动导入分子数据成功。")  
  
  
# ==========================================  
# 布局 2：分析核心展示（仅当数据源就绪后进入）  
# ==========================================  
if df_explore is not None:  
    # 提取表中所有的数值列  
    num_cols = df_explore.select_dtypes(include=[np.number]).columns.tolist()  
    # 剔除可能存在的标识型无用列  
    num_cols = [c for c in num_cols if c not in ["label", "UID", "id", "Unnamed: 0"]]  
    
    if len(num_cols) == 0:  
        st.error("数据校验未通过：未在数据表内扫描出任何数值型物理性质特征字段。")  
    else:  
        # 2.1 配置分析描述符与冗余度自检测  
        st.markdown('<div class="section-header">分析描述符配置与特征冗余诊断</div>', unsafe_allow_html=True)  
        
        with st.container(border=True):  
            col_sel, col_diag = st.columns([0.45, 0.55])  
            
            with col_sel:  
                st.markdown("**选择待考察的自变量描述符**")  
                # 设置默认选项以避免初始空白  
                default_picks = [c for c in ["MW", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds"] if c in num_cols]  
                if not default_picks:  
                    default_picks = num_cols[:min(3, len(num_cols))]  
                    
                selected_cols = st.multiselect(  
                    "请从表格解析出来的字段中，选定需要参与诊断的物性字段：",  
                    options=num_cols,  
                    default=default_picks  
                )  
                
            with col_diag:  
                st.markdown("**相关性共线性检测报告**")  
                if len(selected_cols) > 1:  
                    # 计算特征相关系数矩阵  
                    corr_matrix = df_explore[selected_cols].corr().abs()  
                    redundant_pairs = []  
                    
                    for i in range(len(corr_matrix.columns)):  
                        for j in range(i + 1, len(corr_matrix.columns)):  
                            if corr_matrix.iloc[i, j] > 0.75:  
                                redundant_pairs.append(  
                                    (corr_matrix.columns[i], corr_matrix.columns[j], corr_matrix.iloc[i, j])  
                                )  
                    
                    if redundant_pairs:  
                        warning_html = "<ul style='margin:0; padding-left:20px; font-size:13px; color:#9a3412;'>"  
                        for pair in redundant_pairs:  
                            warning_html += f"<li>特征 <b>{pair[0]}</b> 与 <b>{pair[1]}</b> 相关度过高 (r = {pair[2]:.2f})</li>"  
                        warning_html += "</ul>"  
                        
                        st.markdown(f"""  
                        <div style="background-color: #fffbeb; border-left: 4px solid #d97706; padding: 12px; border-radius: 6px; font-size: 13px; color: #b45309;">  
                            <strong>检测到高共线性冗余特征对：</strong><br>  
                            {warning_html}  
                            <span style="font-size:12px; color:#78350f; display:block; margin-top:8px;">💡 优化建议：共线性特征在机器学习模型训练中会导致权重分配异常，建议在训练分类器前只剔除并只保留高关联组中的其中一个。</span>  
                        </div>  
                        """, unsafe_allow_html=True)  
                    else:  
                        st.markdown("""  
                        <div style="background-color: #f0fdf4; border-left: 4px solid #16a34a; padding: 12px; border-radius: 6px; font-size: 13px; color: #15803d;">  
                            <strong>系统性质自检通过：</strong><br>  
                            当前选定物理特征间的 Pearson 相关性绝对值全数保持在 0.75 以下。物理特征信息隔离良好，可以直接组合用于下游机器学习建模。  
                        </div>  
                        """, unsafe_allow_html=True)  
                else:  
                    st.info("提示：请在左侧多选栏内至少选择 2 个物理描述符，方可激活相关性冗余检测分析。")  
  
  
        # 2.2 特征相关性 & 分布直方图  
        if len(selected_cols) > 0:  
            col_plot_l, col_plot_r = st.columns([0.45, 0.55])  
            
            with col_plot_l:  
                st.markdown('<div class="section-header">特征间相关性系数热图</div>', unsafe_allow_html=True)  
                with st.container(border=True):  
                    if len(selected_cols) >= 2:  
                        fig_corr, ax_corr = plt.subplots(figsize=(5, 4.2))  
                        # 使用 seaborn 渲染干净典雅的蓝红热图  
                        sns.heatmap(  
                            df_explore[selected_cols].corr(),  
                            annot=True,  
                            cmap="RdBu_r",  
                            vmin=-1.0,  
                            vmax=1.0,  
                            fmt=".2f",  
                            annot_kws={"size": 9},  
                            ax=ax_corr  
                        )  
                        ax_corr.tick_params(axis='both', which='major', labelsize=9)  
                        ax_corr.set_title("Pearson Correlation Coefficients", fontsize=10, pad=10)  
                        plt.tight_layout()  
                        
                        st.pyplot(fig_corr)  
                        plt.close(fig_corr)  
                    else:  
                        st.info("相关性热图绘制需要您在配置项中至少挑选 2 个特征。")  
                        
            with col_plot_r:  
                st.markdown('<div class="section-header">物理属性直方分布图</div>', unsafe_allow_html=True)  
                with st.container(border=True):  
                    # 分网格布局平铺直方图，避免图表大小不对称  
                    num_plots = len(selected_cols)  
                    rows = int(np.ceil(num_plots / 2))  
                    cols = 2 if num_plots > 1 else 1  
                    
                    fig_dist, axes = plt.subplots(rows, cols, figsize=(7.5, 2.1 * rows))  
                    
                    # 强行打平成一维数组，方便循环遍历配置  
                    if num_plots == 1:  
                        axes = np.array([axes])  
                    else:  
                        axes = axes.flatten()  
                        
                    for idx, property_name in enumerate(selected_cols):  
                        sns.histplot(  
                            df_explore[property_name].dropna(),  
                            kde=True,  
                            color="#1e3a8a",  
                            ax=axes[idx],  
                            bins=15  
                        )  
                        axes[idx].set_xlabel(property_name, fontsize=8)  
                        
                        # --- 已完成：修改为英文 "Count" ---
                        axes[idx].set_ylabel("Count", fontsize=8)
                        
                        axes[idx].tick_params(labelsize=8)  
                        axes[idx].spines['top'].set_visible(False)  
                        axes[idx].spines['right'].set_visible(False)  
                        axes[idx].set_title(f"{property_name} Density", fontsize=9)  
                        
                    # 如果选择的描述符个数在网格中留空，删除不需要的子图区域  
                    for empty_idx in range(num_plots, len(axes)):  
                        fig_dist.delaxes(axes[empty_idx])  
                        
                    plt.tight_layout()  
                    st.pyplot(fig_dist)  
                    plt.close(fig_dist)  
        else:  
            st.warning("暂无任何特征进入评估，请通过上方的配置栏向清单中添加对应的物理自变量特征。")

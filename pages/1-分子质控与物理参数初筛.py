import streamlit as st   
import pandas as pd   
import numpy as np   
import matplotlib.pyplot as plt   
from rdkit import Chem   
from rdkit.Chem import Descriptors, Lipinski   
import io  
import sys  
import platform

# 1. 解决标准输出流可能存在的编码问题  
if sys.stdout.encoding != 'utf-8':  
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')  

# 2. 更加全面的字体配置，兼容 Windows, macOS, Linux (如 Streamlit Cloud) 
def set_matplot_zh_font():
    # 自动识别系统并配置字体优先级
    plt.rcParams['font.sans-serif'] = [
        'SimHei',             # Windows 常用黑体
        'Microsoft YaHei',    # Windows 微软雅黑
        'Heiti TC',           # macOS 繁体黑体
        'Arial Unicode MS',   # macOS 兼容中文
        'WenQuanYi Micro Hei',# Linux/Ubuntu 下的开源中文字体（部署Streamlit关键）
        'DejaVu Sans',        # 备用英文
        'sans-serif'
    ]
    plt.rcParams['axes.unicode_minus'] = False # 正常显示负号

set_matplot_zh_font()
   
# 设置页面宽屏模式
st.set_page_config(page_title="数据清洗与质控", layout="wide")   
   
# 样式 CSS 保持并微调以适应全宽布局
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

# ==========================================
# 布局重构 1：将过滤配置与文件上传水平放置，节省头部空间
# ==========================================
col_config, col_upload = st.columns([0.45, 0.55])

with col_config:
    st.markdown('<div class="card">', unsafe_allow_html=True)   
    st.subheader("⚙️ 过滤条件配置")   
    f_lip = st.checkbox("必须符合 Lipinski 五规则（允许 1 项超标）", value=True)   
    f_mw = st.checkbox("限制分子量 (MW) 在 160 ~ 500 Da 之间", value=True)   
    f_logp = st.checkbox("限制脂水分配系数 (LogP) 在 -2 ~ 5 之间", value=True)   
    st.markdown('</div>', unsafe_allow_html=True) 

with col_upload:
    st.markdown('<div class="card">', unsafe_allow_html=True)   
    st.subheader("📂 导入分子数据集")   
    csv_file = st.file_uploader("导入分子数据集表 (CSV格式)", type=["csv"], help="表格必须含 smiles 与分类标签 label 列")   
    st.markdown('</div>', unsafe_allow_html=True)   

# ==========================================
# 布局重构 2：文件输入后的预览与核心逻辑全部展开为 “上下结构”
# ==========================================
if csv_file:
    df = pd.read_csv(csv_file)
    
    # --- 1号块: 读入预览块 (占据全宽，不挤在右侧) ---
    st.markdown('<div class="card">', unsafe_allow_html=True)   
    st.subheader("🔍 原始数据集预览")   
    st.info(f"读取到分子共计：{len(df)} 种")   
    st.dataframe(df.head(6), use_container_width=True)   
    st.markdown('</div>', unsafe_allow_html=True)   
    
    if "smiles" not in df.columns:   
        st.error("表格结构异常：缺失了名为 smiles 的分子式列。")   
    else:   
        # 初始化 session_state，持久化数据
        if "cleaned_df" not in st.session_state:  
            st.session_state["cleaned_df"] = None  
            st.session_state["raw_len"] = 0  
            st.session_state["orig_mws"] = []  
            st.session_state["final_mws"] = []  

        # 清洗动作执行排在正中间 (大按钮)
        if st.button("🚀 开始清洗数据", type="primary", use_container_width=True):   
            progress_bar = st.progress(0)   
            status_text = st.empty()   
               
            cleaned_indices = []   
            orig_mws, final_mws = [], []   
            orig_logps, final_logps = [], []   
            raw_len = len(df)   
               
            for idx, row in df.iterrows():   
                if raw_len > 20 and idx % max(1, raw_len // 20) == 0:   
                    val = (idx + 1) / raw_len   
                    progress_bar.progress(min(val, 1.0))   
                    status_text.text(f"质控计算中：{int(val*100)}% ({idx+1}/{raw_len})")   
                   
                s = row["smiles"]   
                try:   
                    # 安全检查
                    if pd.isna(s) or not isinstance(s, str):  
                        continue  
                    s = s.strip()  
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
                        if lip_violations > 1:   
                            continue   
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
            
            # 存入 session_state
            st.session_state["cleaned_df"] = cleaned_df   
            st.session_state["raw_len"] = raw_len  
            st.session_state["orig_mws"] = orig_mws  
            st.session_state["final_mws"] = final_mws  

        # --- 2号块: 结果上下排版区 (清洗完毕后展示) ---
        if st.session_state["cleaned_df"] is not None:  
            cleaned_df = st.session_state["cleaned_df"]  
            raw_len = st.session_state["raw_len"]  
            orig_mws = st.session_state["orig_mws"]  
            final_mws = st.session_state["final_mws"]  

            # 1. 结果表格展示 (全宽)
            st.markdown('<div class="card">', unsafe_allow_html=True)   
            st.subheader("📋 清洗后所得标准数据集 (前10行记录)")   
            st.dataframe(cleaned_df.head(10), use_container_width=True)   
            st.markdown('</div>', unsafe_allow_html=True)   

            # 2. 统计图表区 (由于宽度大，重新排版柱状图长宽比，避免扭曲)
            st.markdown('<div class="card">', unsafe_allow_html=True)   
            st.subheader("📊 过滤前后数据统计屏")   
            
            # 使用大一号且扁平的比例，适合流式上下阅读
            fig, ax = plt.subplots(figsize=(8, 3))   
            bars = ax.barh(["质控前分子数", "质控后分子数"], [raw_len, len(cleaned_df)], color=["#cbd5e1", "#1e3a8a"], height=0.5)   
            ax.set_xlabel("分子统计数量", fontsize=10)   
            
            # 在条形图上添加数值标签
            for bar in bars:
                width = bar.get_width()
                ax.text(width + (raw_len * 0.01), bar.get_y() + bar.get_height()/2, f'{int(width)}', 
                        va='center', ha='left', fontsize=10, fontweight='bold')
            
            ax.spines['top'].set_visible(False)   
            ax.spines['right'].set_visible(False) 
            ax.spines['left'].set_visible(False) 
            plt.tight_layout()
            
            # 使用 Streamlit 渲染 Matplotlib 图片
            st.pyplot(fig)   
            plt.close(fig) # 释放内存
            st.markdown('</div>', unsafe_allow_html=True)   
               
            # 3. 诊断报告 (全宽)
            avg_mw_diff = np.mean(final_mws) - np.mean(orig_mws) if final_mws and orig_mws else 0   
            st.markdown(f"""   
            <div class="report-card">   
                <h4 style="margin-top:0px; color:#1e3a8a !important;">🔬 数据清洗质控诊断报告</h4>   
                <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">   
                    数据质控流程判定完毕。系统已从原数据库的 {raw_len} 个分子中清除了 <strong>{raw_len - len(cleaned_df)}</strong> 个不合规分子（数据通过保留率：<strong>{(len(cleaned_df)/raw_len * 100):.2f}%</strong>）。   
                    质控后整体小分子均分子量发生 <strong>{avg_mw_diff:.2f} Da</strong> 的位移。<br>   
                    移去物理化学常数过大、或已知存在结构错误的配体原子体系，对于避免后续的建模过拟合具备重要保障意义。   
                </p>   
            </div>   
            """, unsafe_allow_html=True)   
               
            # 4. 下载按钮 (底部，占据主屏)
            csv = cleaned_df.to_csv(index=False).encode('utf-8-sig')   
            st.markdown("<br>", unsafe_allow_html=True)   
            st.download_button("💾 导出质控后的分子数据集 (CSV)", csv, "QC_Cleaned_Dataset.csv", use_container_width=True)   
else:   
    st.info("💡 提示：请先在栏目上方上传含有 SMILES 结构的 CSV 数据表。")

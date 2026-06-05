import streamlit as st
import pandas as pd
import numpy as np
import joblib
import io
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score,
                             roc_auc_score, roc_curve)

# 字体配置，确保中文不乱码
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="QSAR预测建模", layout="wide")

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
    .stButton>button {
        background-color: #1e3a8a !important;
        color: white !important;
        border-radius: 8px !important;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>QSAR 分子活性分类预测建模</h1>
    <p>提供多重主流算法的对比学习，拟合小分子的化学结构与其活性关系，并提取对活性贡献最高的物理特征参数。</p>
</div>
""", unsafe_allow_html=True)

# 独立的数据源选择器
st.markdown('<div class="card">', unsafe_allow_html=True)
st.subheader("数据源选择")
data_source = st.radio(
    "请指定建模训练所解析的数据表：",
    ["使用前序质控清洗的文件缓存", "上传外部包含活性标记 (label) 的 CSV 表格"]
)

df_combined = None
if data_source == "使用前序质控清洗的文件缓存":
    if "cleaned_df" in st.session_state:
        df_combined = st.session_state["cleaned_df"]
        st.success("成功载入前序清洗过的缓存数据。")
    else:
        st.error("缓存中不存在数据。请先执行过滤清洗或选择上传包含活性 label 的本地分类 CSV 数据集。")
else:
    up_csvs = st.file_uploader("导入本地模型训练集 (表内须包含 label 列以区分类别)", type=["csv"])
    if up_csvs:
        df_combined = pd.read_csv(up_csvs)
st.markdown('</div>', unsafe_allow_html=True)

if df_combined is not None:
    preset_phys = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds"]
    valid_features = [c for c in preset_phys if c in df_combined.columns]
    
    morgan_cols = [col for col in df_combined.columns if col.startswith("MorganFP")]
    all_features = valid_features + morgan_cols

    if len(all_features) == 0 or "label" not in df_combined.columns:
        st.error("字段错误：上传的数据表中必须包含 label 列作为活性分类标记（1 为有活性，0 为无活性），并且要有物理指标参数。")
    else:
        col_ctrl, col_view = st.columns([0.35, 0.65])
        
        with col_ctrl:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("建模参数设置")
            ratio = st.slider("测试集切分划分比例 (%)", min_value=10, max_value=40, value=20)
            btn_calc = st.button("启动对比训练", type="primary", use_container_width=True)
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("物理特征矩阵属性")
            st.write(f"- 自变量参与建模维度：{len(all_features)}")
            st.write(f"- 具有生物活性的实验样本数 (label=1)：{df_combined['label'].sum()} / {len(df_combined)}")
            st.markdown('</div>', unsafe_allow_html=True)
            
        with col_view:
            if btn_calc:
                X = df_combined[all_features].fillna(df_combined[all_features].mean())
                y = df_combined["label"]
                
                scaler = StandardScaler()
                X_scaled = scaler.fit_transform(X)
                
                X_train, X_test, y_train, y_test = train_test_split(
                    X_scaled, y, test_size=ratio/100.0, random_state=42, stratify=y
                )
                
                models = {
                    "随机森林 (Random Forest)": RandomForestClassifier(n_estimators=100, max_depth=8, random_state=42),
                    "支持向量机 (Support Vector Machine)": SVC(C=1.0, probability=True, random_state=42),
                    "逻辑回归 (Logistic Regression)": LogisticRegression(C=1.0, max_iter=300, random_state=42)
                }
                
                perf_list = []
                # 绘制中文 ROC 曲线图
                fig_roc, ax_roc = plt.subplots(figsize=(4.8, 3.8))
                trained_forest = None
                
                for name, model in models.items():
                    model.fit(X_train, y_train)
                    y_pred = model.predict(X_test)
                    y_pred_proba = model.predict_proba(X_test)[:, 1]
                    
                    acc = accuracy_score(y_test, y_pred)
                    prec = precision_score(y_test, y_pred)
                    rec = recall_score(y_test, y_pred)
                    f1 = f1_score(y_test, y_pred)
                    auc_val = roc_auc_score(y_test, y_pred_proba)
                    
                    perf_list.append({
                        "模型算法": name,
                        "准确率": round(acc, 3),
                        "精确度": round(prec, 3),
                        "召回率": round(rec, 3),
                        "F1评分": round(f1, 3),
                        "AUC分类面积": round(auc_val, 3)
                    })
                    
                    # 避免 matplotlib 中文报错，线标签可使用英文
                    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
                    ax_roc.plot(fpr, tpr, label=f"{name.split(' (')[0]} (AUC={auc_val:.3f})")
                    
                    if "Random Forest" in name or "随机森林" in name:
                        trained_forest = model
                        
                ax_roc.plot([0, 1], [0, 1], color="#94a3b8", linestyle="--")
                ax_roc.set_xlabel("假阳性率 (FPR)", fontsize=8)
                ax_roc.set_ylabel("真阳性率 (TPR)", fontsize=8)
                ax_roc.set_title("各算法活性预测 ROC 曲线对标", fontsize=9)
                ax_roc.legend(loc="lower right", fontsize=8)
                ax_roc.spines['top'].set_visible(False)
                ax_roc.spines['right'].set_visible(False)
                
                st.markdown('<div class="card">', unsafe_allow_html=True)
                col_p_l, col_p_r = st.columns([0.55, 0.45])
                with col_p_l:
                    st.subheader("机器学习模型 ROC 曲线")
                    st.pyplot(fig_roc)
                with col_p_r:
                    st.subheader("不同算法精度评比榜单")
                    df_perf = pd.DataFrame(perf_list).sort_values("AUC分类面积", ascending=False)
                    st.dataframe(df_perf, index=False, use_container_width=True)
                    
                    best_name = df_perf.iloc[0]["模型算法"]
                    buf = io.BytesIO()
                    joblib.dump(models[best_name], buf)
                    buf.seek(0)
                    st.markdown("<br>", unsafe_allow_html=True)
                    st.download_button(
                        label=f"将最佳模型 ({best_name.split(' (')[0]}) 导出 (.pkl)",
                        data=buf,
                        file_name="qsar_best_model.pkl",
                        mime="application/octet-stream",
                        use_container_width=True
                    )
                st.markdown('</div>', unsafe_allow_html=True)
                
                # 特征重要性卡片
                st.markdown('<div class="card">', unsafe_allow_html=True)
                st.subheader("描述符特征解释贡献分析")
                ch_map = {
                    "MW": "分子量 (MW)", 
                    "LogP": "脂水分配系数 (LogP)", 
                    "TPSA": "极性拓扑表面积 (TPSA)",
                    "HBD": "氢键供体数 (HBD)", 
                    "HBA": "氢键受体数 (HBA)", 
                    "RotatableBonds": "柔性单键数 (RB)"
                }
                
                if trained_forest is not None:
                    importances = trained_forest.feature_importances_
                    indices = np.argsort(importances)[::-1][:12]
                    top_features = [all_features[i] for i in indices]
                    top_importances = importances[indices]
                    
                    c_bar_l, c_bar_r = st.columns([0.55, 0.45])
                    with c_bar_l:
                        fig_bar, ax_bar = plt.subplots(figsize=(5.5, 3.8))
                        # 使用特征中文名绘图，排除中文乱码
                        features_zh_labels = [ch_map.get(x, x) for x in top_features]
                        ax_bar.barh(features_zh_labels[::-1], top_importances[::-1], color="#1e3a8a", edgecolor="#0f172a", height=0.55)
                        ax_bar.set_xlabel("自变量相对贡献度 (Gini重要性)", fontsize=9)
                        ax_bar.set_title("调控活性关键理化属性排名", fontsize=9)
                        ax_bar.spines['top'].set_visible(False)
                        ax_bar.spines['right'].set_visible(False)
                        st.pyplot(fig_bar)
                    with c_bar_r:
                        primary_feat = top_features[0]
                        primary_feat_zh = ch_map.get(primary_feat, primary_feat)
                        st.markdown(f"""
                        <div class="report-card" style="margin-top:0px;">
                            <h4 style="margin-top:0px; color:#1e3a8a !important;">QSAR 构效关系分析报告</h4>
                            <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                                经树结构信息增益判定：<strong>{primary_feat_zh}</strong> 是调控当前活性类别最核心的分子物理参数。其贡献度达到了 <strong>{(top_importances[0]*100):.1f}%</strong>。<br><br>
                                <strong>结构微调方向：</strong><br>
                                - 建议在继续设计子代衍生物时，密切监控该特征的区间变动。如果性质是 <strong>LogP</strong>，表明靶点结合腔具有很强的不带电疏水匹配倾向，调节分子周围脂肪性碳链的支链排布可有效控制化合物的结合活性。
                            </p>
                        </div>
                        """, unsafe_allow_html=True)
                st.markdown('</div>', unsafe_allow_html=True)
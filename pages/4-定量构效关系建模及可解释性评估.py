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

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

st.set_page_config(page_title="QSAR模型训练", layout="wide")

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
        margin-top: 15px;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="blue-banner">
    <h1>QSAR 生物活性分类预测模型训练</h1>
    <p>比对不同主流机器学习算法的拟合性能，绘制分类 ROC 评估曲线，提取影响分子活性的关键理化描述符。</p>
</div>
""", unsafe_allow_html=True)

st.markdown('<div class="card">', unsafe_allow_html=True)
st.subheader("数据源选择")
data_source = st.radio(
    "指定建模训练所解析的数据表：",
    ["调取前序质控清洗的文件缓存", "上传外部包含活性标记 (label) 的分类 CSV 库表"]
)

df_target = None
if data_source == "调取前序质控清洗的文件缓存":
    if "cleaned_df" in st.session_state:
        df_target = st.session_state["cleaned_df"]
        st.success("读取成功：已载入过滤器缓存分子。")
    else:
        st.error("系统没有检测到清洗数据缓存。请先进行第一步清洗过滤，或者在上方修改为上传本地 CSV。")
else:
    f_csv = st.file_uploader("上传本地训练集 CSV (表格须含有 label 分类列且为 0 或 1)", type=["csv"])
    if f_csv:
        df_target = pd.read_csv(f_csv)
st.markdown('</div>', unsafe_allow_html=True)

if df_target is not None:
    preset_phys = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds"]
    trainable_features = [col for col in preset_phys if col in df_target.columns]
    
    if len(trainable_features) == 0 or "label" not in df_target.columns:
        st.error("字段异常：建模特征集或目标活性标记列 'label' 缺失！")
    else:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("模型配置版")
        col_t1, col_t2 = st.columns(2)
        with col_t1:
            ratio = st.slider("测试集切分划分比例 (%)", min_value=15, max_value=40, value=25)
        with col_t2:
            st.write(f"- 自变量参与建模维度数：{len(trainable_features)}")
            st.write(f"- 具有生物活性样本分子数 (label=1)：{df_target['label'].sum()} / {len(df_target)}")
            
        btn_train = st.button("启动本轮建模对比训练", type="primary", use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)
        
        if btn_train:
            X = df_target[trainable_features].fillna(df_target[trainable_features].mean())
            y = df_target["label"]
            
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            
            X_train, X_test, y_train, y_test = train_test_split(
                X_scaled, y, test_size=ratio/100.0, random_state=42, stratify=y
            )
            
            models = {
                "随机森林 (Random Forest)": RandomForestClassifier(n_estimators=100, max_depth=7, random_state=42),
                "支持向量机 (Support Vector Machine)": SVC(C=1.0, probability=True, random_state=42),
                "逻辑回归 (Logistic Regression)": LogisticRegression(C=1.0, max_iter=300, random_state=42)
            }
            
            perf_list = []
            fig_roc, ax_roc = plt.subplots(figsize=(6, 4.2))
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
                    "F1得分": round(f1, 3),
                    "AUC面积": round(auc_val, 3)
                })
                
                fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
                ax_roc.plot(fpr, tpr, label=f"{name.split(' (')[0]} (AUC={auc_val:.3f})")
                
                if "Random Forest" in name:
                    trained_forest = model
                    
            ax_roc.plot([0, 1], [0, 1], color="#94a3b8", linestyle="--")
            ax_roc.set_xlabel("FPR (假阳性率)", fontsize=9)
            ax_roc.set_ylabel("TPR (真阳性率)", fontsize=9)
            ax_roc.set_title("各分类器测试集 ROC 对标曲线", fontsize=10)
            ax_roc.legend(loc="lower right", fontsize=8)
            ax_roc.spines['top'].set_visible(False)
            ax_roc.spines['right'].set_visible(False)
            
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("模型计算拟合性能评估")
            
            col_l, col_r = st.columns([0.5, 0.5])
            with col_l:
                st.pyplot(fig_roc)
            with col_r:
                df_perf = pd.DataFrame(perf_list).sort_values("AUC面积", ascending=False)
                st.dataframe(df_perf, use_container_width=True)
                
                best_model_name = df_perf.iloc[0]["模型算法"]
                best_model = models[best_model_name]
                
                buf = io.BytesIO()
                joblib.dump(best_model, buf)
                buf.seek(0)
                st.write("")
                st.download_button(
                    label=f"直接下载表现最佳的分类模型 ({best_model_name.split(' (')[0]}) (.pkl)",
                    data=buf,
                    file_name="qsar_classifier_model.pkl",
                    mime="application/octet-stream",
                    use_container_width=True
                )
            st.markdown('</div>', unsafe_allow_html=True)
            
            if trained_forest is not None:
                st.markdown('<div class="card">', unsafe_allow_html=True)
                st.subheader("模型物理化学描述符自变量贡献度")
                
                ch_map = {
                    "MW": "分子量 (MW)", 
                    "LogP": "水油分配系数 (LogP)", 
                    "TPSA": "极性拓扑面积 (TPSA)",
                    "HBD": "氢键供体数 (HBD)", 
                    "HBA": "氢键受体数 (HBA)", 
                    "RotatableBonds": "可旋转单键 (RB)"
                }
                
                importances = trained_forest.feature_importances_
                indices = np.argsort(importances)[::-1]
                top_features = [trainable_features[i] for i in indices]
                top_importances = importances[indices]
                
                col_bar, col_card = st.columns([0.5, 0.5])
                with col_bar:
                    fig_bar, ax_bar = plt.subplots(figsize=(5.5, 3.4))
                    features_zh = [ch_map.get(f, f) for f in top_features]
                    ax_bar.barh(features_zh[::-1], top_importances[::-1], color="#1e3a8a", edgecolor="#0f172a", height=0.5)
                    ax_bar.set_xlabel("相对重要性百分比 (Gini权重)")
                    ax_bar.set_title("物理属性解释重要度排行")
                    ax_bar.spines['top'].set_visible(False)
                    ax_bar.spines['right'].set_visible(False)
                    st.pyplot(fig_bar)
                with col_card:
                    core_feat = top_features[0]
                    st.markdown(f"""
                    <div class="report-card" style="height: 100%;">
                        <h4 style="margin-top:0px; color:#1e3a8a !important;">QSAR 关联解析结论</h4>
                        <p style="color:#475569; font-size:14px; line-height:1.6; margin:0;">
                            通过随机森林分类算法挖掘得知，<strong>{ch_map.get(core_feat, core_feat)}</strong> 是调控当前这批活性状态最为关键的理化性质，主导了约 <strong>{(top_importances[0]*100):.1f}%</strong> 的解释特征。<br><br>
                            - <strong>构效改良机理指引：</strong> 如果关键决策项为脂溶常数 (LogP)，提示目标受体可能存在较强的疏水结合倾向。调控化合物上的非极性支链长度将能够比较有效地控制其活性分布状态。
                        </p>
                    </div>
                    """, unsafe_allow_html=True)
                st.markdown('</div>', unsafe_allow_html=True)

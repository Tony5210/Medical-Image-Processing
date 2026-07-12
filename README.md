# Spatial Heterogeneity Analysis of White Matter Hyperintensities (WMH)
# 脑白质高信号(WMH)空间异质性分析

**Medical Image Processing Course - Final Project**
**医学图像处理课程 - 期末大作业**

### 目录
main.m --- 源代码\n
MIAreport.pdf --- 报告\n
Fig1_Registration_Results.png\n
Fig2_Clustered_WMH_Regions.png\n
WMH_Regional_Associations.png\n
--- 图片

### 📌 项目简介

脑白质高信号（WMH）是反映脑健康状况及神经退行性疾病（如阿尔茨海默病）的重要影像学表型。然而，不同个体的 WMH 在空间分布上表现出显著的异质性。

本项目实现了一套全自动的神经影像处理流水线，主要任务包括：

1. **图像配准**：将 T1 和 T2 FLAIR 图像配准并空间归一化至 MNI 标准空间。
2. **空间聚类**：利用 K-Means 算法对全脑 WMH 进行降维，划分为具有解剖学意义的优势区域。
3. **关联分析**：计算各区域的病灶占比，并分析其与临床表型数据（如年龄、性别、APOE4、ADAS11）的皮尔逊相关性。

### ⚙️ 方法与技术路线

完整的处理流程封装于单一脚本（`main.m`）中，具体步骤如下：

* **图像配准 (基于 SPM12)**: 
  * 将 T2 FLAIR 刚性配准至原生 T1 空间（使用标准化互信息 NMI）。
  * 将 T1 图像非线性配准至 MNI152 标准模板。
  * 将 WMH 二值掩码空间映射至 MNI 空间（强制使用最近邻插值）。
* **掩码后处理与聚类**: 
  * 三维高斯平滑与空间降采样。
  * 基于体素特征的 K-Means 聚类降维（选取 $k=4$）。
* **表型关联**: 
  * 计算区域病灶占比，并绘制包含线性回归与显著性高亮的散点图。

### 🚀 环境依赖与使用说明

#### 运行要求

* **MATLAB** (推荐 R2023b 及以上版本)
* **SPM12** 工具箱（需添加至 MATLAB 路径）。

#### 运行方法

1. 克隆本仓库。

2. 在 MATLAB 中配置 SPM12 路径：

   ```matlab
   addpath('/path/to/spm12');
   ```

3. 由于 ADNI 数据使用协议限制，本仓库**未包含**原始 NIfTI 影像数据及临床数据表（`all30m.xlsx`）。若需运行代码，请通过 [ADNI 官网](https://adni.loni.usc.edu/) 获取授权，并将数据放置于项目要求的 `images/` 目录下。

4. 在命令窗口运行主脚本：

   ```matlab
   main
   ```

### ⚠️ 免责声明与数据隐私

本项目使用了来自阿尔茨海默病神经影像学计划（ADNI）的数据。根据 ADNI 数据使用协议（DUA）的要求，**本仓库绝不包含、也不公开分享任何个体级别的原始数据、衍生影像或临床表型文件**。本仓库内开源的代码仅用于教育学习、方法学展示及同行交流。

---

### 📌 Project Overview
White Matter Hyperintensities (WMH) are key imaging biomarkers for cerebrovascular and neurodegenerative diseases (e.g., Alzheimer's disease). However, WMH exhibits significant spatial heterogeneity across individuals. 

This project implements a fully automated neuroimaging pipeline to:
1. **Coregister** and **Normalize** T1 and T2 FLAIR images to the MNI standard space.
2. **Cluster** spatial WMH distributions into anatomically meaningful regions using K-Means.
3. **Analyze** the statistical correlation between regional lesion burdens and clinical phenotypic data (e.g., Age, Gender, APOE4, ADAS11).

### ⚙️ Methodology & Pipeline
The entire workflow is encapsulated in a single script (`main.m`) and consists of the following steps:
* **Image Registration (SPM12)**: 
  * Rigid co-registration of T2 FLAIR to native T1 (Normalized Mutual Information cost function).
  * Non-linear normalization of T1 to the MNI152 template.
  * Spatial mapping of WMH binary masks to the MNI space (Nearest Neighbor interpolation).
* **Post-processing & Clustering**: 
  * 3D Gaussian smoothing and spatial downsampling.
  * Voxel-wise dimensionality reduction via K-Means clustering ($k=4$).
* **Phenotypic Association**: 
  * Calculation of regional lesion proportions and Pearson correlation with clinical phenotypes.

### 🚀 Prerequisites and Usage
#### Requirements
* **MATLAB** (Tested on R2023b+)
* **SPM12** (Statistical Parametric Mapping) added to the MATLAB path.

#### How to Run
1. Clone this repository.
2. Ensure SPM12 is installed and initialized in your MATLAB environment:
   ```matlab
   addpath('/path/to/spm12');
   ```
3. Due to ADNI Data Use Agreements, the raw NIfTI images and clinical spreadsheets (`all30m.xlsx`) are **not** included in this repository. To run the code, you must obtain authorized access from the [ADNI Database](https://adni.loni.usc.edu/) and place the data in the `images/` directory as expected by the script.
4. Run the main script:
   ```matlab
   main
   ```

### ⚠️ Disclaimer and Data Privacy
This project uses data from the Alzheimer's Disease Neuroimaging Initiative (ADNI). Per the ADNI Data Use Agreement (DUA), **no participant-level raw data, derived images, or clinical phenotypes are shared in this repository**. The uploaded code is intended solely for educational, methodological demonstration, and peer-review purposes.

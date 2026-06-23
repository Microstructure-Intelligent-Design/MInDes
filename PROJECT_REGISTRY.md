# MInDes 项目认知与开发日志

> **文件用途**：该文件从 AI 协作视角记录 MInDes 项目的完整认知基线、架构决策、开发历史与当前状态。
> **每次新会话开始时应首先阅读此文件**，以快速重建项目上下文，避免重复探索。
> **每次会话结束时，开发者应更新本文件**，将本次会话的认知增量沉淀下来。

---

## 一、项目身份卡

| 字段 | 内容 |
|------|------|
| **项目全称** | Microstructure Intelligent Design Core（MInDes） |
| **定位** | 相场法微结构演化求解器（Phase-Field Simulation Solver） |
| **语言** | C++17（主），Python 3（辅助工具） |
| **许可证** | GPLv3 |
| **开发者联系** | 黄奇博士 — qihuang0908@163.com |
| **仓库根目录** | E:\GITHUB\Microstructure-Intelligent-Design\MInDes |
| **实际源码目录** | MInDes/MInDes/src/ |
| **解决方案文件** | MInDes/MInDes.sln（VS2022） |
| **构建平台** | Windows MSVC v143（主）/ Linux g++ + OpenMP |
| **关键依赖** | FFTW 3，预编译静态库（libMInDeslib.a, libOptAlgorithm.a） |

---

## 二、架构全景

### 2.1 三层模块化流水线

main.cpp 入口 -> MainIterator 主循环调度器

三个阶段：
- **init_modules()**: 解析 .mindes -> 注册模块 -> 创建输出目录
- **run()**:
  - Pre-Exec(i/ii/iii): 预处理（微结构初始化、填充合并）
  - Exec(i/ii/iii): 主循环计算（DDC/DDCCPAI/GGS 物理模型）
  - Pos-Exec(i/ii/iii): 后处理（VTS输出、力学场、流体LBM、ML推理）
  - Deinit: 清理、日志汇总

### 2.2 模块注册机制

Module.h 定义 Solver_Module 结构体（函数指针集合），通过 load_a_new_module() 注册到全局 module_list。每个阶段分 i/ii/iii 三个子阶段控制执行序。

### 2.3 关键源码文件映射

| 层次 | 路径（相对 MInDes/MInDes/src/） | 职责 |
|------|------|------|
| 入口 | main.cpp | 入口 |
| 调度 | MainIterator.h/.cpp | 主循环调度器 |
| 基础 | modules/base/MACRO_DEF.h | 类型别名（REAL=double）、跨平台宏、工具函数 |
| 基础 | modules/base/VectorMatrix.h/.cpp | 向量/矩阵运算 |
| 基础 | modules/base/RotationMatrix.h/.cpp | 旋转矩阵（欧拉角） |
| 基础 | modules/base/Mesh_0.h | 网格数据结构 |
| 输入 | modules/input_modules/inputfiles/InputFileReader.h/.cpp | .mindes 文件解析 |
| 输入 | modules/input_modules/inputfiles/UserStartUp.h/.cpp | 用户交互启动 |
| 预处理 | modules/preprocess_modules/MicrostructureInit.h/.cpp | 微结构初始化 |
| 预处理 | modules/preprocess_modules/MicrostructureInit/VoronoiStructure.h/.cpp | Voronoi 晶粒 |
| 预处理 | modules/preprocess_modules/MicrostructureInit/PorousStructure.h/.cpp | 多孔结构 |
| 预处理 | modules/preprocess_modules/MicrostructureInit/Bmp24Structure.h/.cpp | BMP 图片导入 |
| 预处理 | modules/preprocess_modules/MicrostructureInit/GeometryStructure.h | 几何体定义 |
| 预处理 | modules/preprocess_modules/Pretreatment.h/.cpp | 填充、合并、重排序 |
| 模型 | modules/model_modules/grain_grows_spinodal/ | 晶粒生长+旋节分解（GGS） |
| 模型 | modules/model_modules/Data_Driven_Complex/ | 数据驱动复杂模型（DDC） |
| 模型 | modules/model_modules/ddc_calphad_ai/ | DDC+CALPHAD+AI（DDCCPAI） |
| 后处理 | modules/postprocess_modules/WriteVTS.h/.cpp | VTK/VTS 输出 |
| 后处理 | modules/postprocess_modules/MechanicalField.h/.cpp + Mechanics/ | 力学场（弹性/塑性） |
| 后处理 | modules/postprocess_modules/FluidField.h/.cpp + FluidDynamics/ | 流体 LBM |
| 后处理 | modules/postprocess_modules/MachineLearning.h/.cpp + FieldPredictor/ | ML 推理（共享内存 IPC） |
| 后处理 | modules/postprocess_modules/AutoDeltTime.h/.cpp | 自适应时间步长 |
| 工具 | tools/ML_field_predictor/service.py | Python ML 推理服务 |
| 工具 | tools/ML_field_predictor/shm_layout.py | 共享内存布局定义 |

---

## 三、认知基线

### 3.1 核心领域概念

- **相场法**：通过连续序参量描述微结构演化，自由能泛函变分推导 Allen-Cahn / Cahn-Hilliard 方程
- **CALPHAD**：基于热力学数据库计算相图与热力学性质
- **DDC**: Data-Driven Complex，数据驱动方式构造演化方程
- **DDCCPAI**: DDC + CALPHAD + AI 三层耦合
- **GGS**: Grain Growth & Spinodal decomposition
- **LBM**: Lattice Boltzmann Method 流体动力学

### 3.2 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 数值类型 | REAL = double | 科学计算精度 |
| 并行方案 | OpenMP | 共享内存多线程 |
| 网格结构 | 结构化正交网格 | 有限差分法求解 PDE |
| 输入格式 | 自定义 .mindes（类似 INI） | 灵活键值对+数学表达式 |
| 模块注册 | 函数指针向量 | 轻量零抽象开销 |
| 跨平台 | #ifdef 条件编译 | 无额外抽象层 |
| ML推理IPC | Windows共享内存（CreateFileMapping） | 低延迟 |

### 3.3 编码约定

- 命名空间: pf，子命名空间分层（如 pf::main_iterator）
- 文件名: PascalCase
- 文件头: #pragma once
- 类型别名: REAL = double

---

## 四、会话日志

按时间倒序排列。

---

### 会话 #002 — 2026-06-22

**目标**：扩展 DDCCPAI 模块，在 Pre-Exec III 阶段生成可配置的能量最小化 CSV 扫描数据。

**已完成的工作**：
- 实现 `write_csv_energy_minimization_results()`，支持组分、region 内全部相分数和温度的规则扫描
- 增加 `.mindes` CSV 扫描参数解析及 region、phase、concentration、范围和步长校验
- 限制每个组分位于 `[0,1)` 且样本总组分严格小于 1；每个相分数大于 `PhiCon_Cut_Off` 且 region 相分数和严格小于 1
- CSV 输出总组分、相分数、平衡相组分、化学势、化学能和收敛信息，并逐样本输出详细进度日志
- 增加独立的 `Model.DDCCPAI.Output.CSV.ChemEnergy` 相选择、组分/温度扫描参数，为每个相独立输出 `parameters::fchem` 化学能 CSV
- 在 DDCCPAI `exec_pre_iii()` 中无条件调用 CSV 函数，启用判断封装在函数内部
- MSVC x64/Release 构建验证通过

---

### 会话 #001 — 2026-06-22

**目标**：首次接触项目，全量文件扫描与认知梳理，创建日志文件。

**已完成的工作**：
- 遍历扫描全部项目文件（约 50+ 编译单元，1.5 万+ 行源码）
- 深入阅读核心源码：main.cpp -> MainIterator -> Module.h -> MACRO_DEF.h -> 各模块接口
- 理解三层模块化流水线架构
- 理解 .mindes 输入文件格式与解析流程
- 识别所有物理模型模块（GGS / DDC / DDCCPAI / LBM / Mechanics）
- 识别 ML 推理子系统的共享内存 IPC 架构
- 创建本 PROJECT_REGISTRY.md 文件

**关键发现**：
- 代码库约 50+ 个编译单元
- x64/Debug/ 和 x64/Release/ 存在大量 .obj 中间产物
- 部分文件仅为头文件实现（Mesh_0.h, GeometryStructure.h, DDC_Manager.h 等）
- Linux 构建依赖 libMInDeslib.a、libOptAlgorithm.a 两个预编译库（未纳入本仓库）
- .gitignore 排除 MInDes/MInDes/infile/ 但保留示例目录

**决策**：
- 创建 PROJECT_REGISTRY.md 用于后续会话的上下文重建

---

## 五、当前状态

| 维度 | 状态 |
|------|------|
| **代码可编译** | 存在 Debug/Release 构建产物，推测可编译 |
| **代码完整性** | 全部源文件已梳理，各模块接口已理解 |
| **文档完整性** | 有 README + CONTRIBUTING，无 API/设计文档 |
| **测试** | 无单元测试或集成测试框架 |
| **活跃开发** | 未知 |

---

## 六、已知技术债务

| 编号 | 问题 | 严重度 |
|------|------|--------|
| TD-001 | OpenMP Debug 模式下强制单线程 | 低 |
| TD-002 | whereami.c 是 C 文件混入 C++ 项目 | 低 |
| TD-003 | .obj 中间文件本地磁盘占用大 | 低 |
| TD-004 | libMInDeslib.a / libOptAlgorithm.a 为闭源预编译库 | 中 |
| TD-005 | 跨平台路径依赖 dirSeparator 宏，拼接处可能遗漏 | 低 |

---

## 七、术语表

| 缩写 | 全称 | 含义 |
|------|------|------|
| phi | Phase-field variable | 相场序参量 |
| Con (c) | Concentration | 浓度场 |
| DDC | Data-Driven Complex | 数据驱动复杂模型 |
| DDCCPAI | DDC+CALPHAD+AI | 三层耦合模型 |
| GGS | Grain Growth & Spinodal | 晶粒生长与旋节分解 |
| LBM | Lattice Boltzmann Method | 格子玻尔兹曼法 |
| VTS | VTK Structured Grid | 可视化格式 |
| CALPHAD | CALculation of PHAse Diagrams | 相图计算 |

---

## 八、附录

### A. 构建指引

**Windows (MSVC 2022)**: 打开 MInDes/MInDes.sln -> x64/Release -> 生成
**Linux (g++)**: cd MInDes/MInDes && bash scripts/compile_exec.sh

### B. .mindes 语法要点

```
Solver.Mesh.Nx = 64
Solver.Mesh.PCT = (phi_number, con_number, is_temperature_on)
Preprocess.Microstructure.geometry_layer_0.property = (phi_index, geo_type, rot_gauge, reverse)
Solver.Output.VTS.frequence = 1000
```

---

> **维护者**：每次会话的 AI 或开发者协作更新。
> **版本**：2026-06-22 v1

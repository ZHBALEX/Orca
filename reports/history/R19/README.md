# Orca 涡命名核查：方法与文献

> 历史报告；以最新总报告为准。原报告散落附件不随本合集上传。

历史报告快照；部分判断后续已更新。请先阅读最新总报告。原始数据附件保留本地，本合集不包含数据或论文全文。

返回精简结论 [original local attachment: Orca_涡命名与新颖性核查.html]

# 方法、证据边界与文献核查

## 原始数据与显示口径

分析 baseline / Sweep35、isolated / Sweep35_OnlyFin、body_only / NoFin。三例各 17 帧，共 51 帧。baseline 与 NoFin 步数 4320–4800；isolated 为 1440–1920；每 30 步一帧，周期 480 步。每例 q-fin-top 原图为 16 张，从相位 0.0625 到 1；当前标作相位 0 的图通常用周期末帧代替，数值分析同时保留首尾帧。三个算例端点 Q 投影相对差分别为 baseline: 0.00179, isolated: 0.0585, body_only: 0.000551，不把全周期所有统计重复计入端点。

局部体积 X=[10.10,10.63]，Y=[5.70,6.06]，Z=[5.72,5.98]，针对负 Z 侧胸鳍和近尾迹。速度梯度完全从流体 stencil 计算；固体附近额外排除两层。以 Q=−½ tr[(∇u)²] 重算并与 Q/Wz 比较，51 帧最低相关 0.999999961、最大相对 L2 误差 0.000266。上游平均 U 均为正：0.2801–0.2916。Q=50/250 来自现存样式；没有把原图蓝灰色当成旋向。

## 身份追踪与不确定性

分别在 Z/L=5.82、5.86、5.89 截取原网格，用 Q 阈值 25、50、100、250 和带符号 ωz 分离旋转核心；少于 5 个截面单元的结构不作跟踪。相邻帧按局部速度预测位置，以同符号最小距离进行一对一匹配，距离门槛 0.015L 或 0.023L；不跨越缺失帧。该方法建立 Eulerian 特征对应，不能证明材料涡管身份、合并或重联。

负 ωz 核在 baseline 和 isolated 的 0.5625–1.0 连续输运，可在三处截面重复看到。中段截面 Q=25、50 下对两种门槛均稳健；高阈值可能切断轨迹，不能把轨迹断裂当作涡生成／消失。部分外展向轨迹也持续到 Q=100。正负尾迹均存在；单张图中的 TEV₂ 不能仅因位于尾缘下游就定来源。3D Q 连通分量另外保留，但阈值把同一物理涡拆开或把不同涡连在一起，因此不按连通分量数量宣称新涡数。

代表截面核心追踪 [original local attachment: section_core_correspondence.json] · 三截面、四阈值、两门槛结果 [original local attachment: track_sensitivity.json] · 原图清单与显示核查 [original local attachment: display_and_data_audit.json]

## 鳍根区定量比较

固定区域 X=[10.32,10.43]、Y=[5.89,5.975]、Z=[5.89,5.945]；每个算例在自身流体单元上积分，使用局部网格间距估算单元体积。表中为 16 个不重复相位的均值 ∫Q·1(Q>Q₀)dV。该区含多个涡／剪切分支，而且身体会改变可用流体体积，不能当作同一材料涡的强度比或力份额。

| Q₀ | baseline | isolated | No Fin | baseline / isolated |
| --- | --- | --- | --- | --- |
| 25 | 0.00196435 | 0.000856065 | 8.32526e-05 | 2.29 |
| 50 | 0.00160615 | 0.000394341 | 2.81237e-05 | 4.07 |
| 100 | 0.00103348 | 3.15543e-05 | 7.91205e-06 | 32.75 |
| 250 | 0.000239193 | 0 | 0 | 不定义 |

三个非零对照阈值下均有增强，但倍数高度阈值敏感；Q=250 时 isolated 区域中没有保留单元，不能报告无限倍增强。baseline 的强核与反旋分支是可重复观测，独立新生涡及其壁面来源仍没有被确认。

照 Liu et al. (2017) 尝试以峰值涡量的 10% 闭合线计算截面环量，并另查 20%、40%。X≈10.35、t/T=0 的 baseline 正核峰值 ωx=34.25；10% 区域与排除带相接，不闭合；20% 线闭合，Γx=0.006805。isolated 同一搜索框的最强正核位置不同，不能把两者 Γ 比值当成同一涡的增强。t/T=0.5 的 baseline 10%–40% 轮廓均接触排除带，不能提供可靠独立涡环量。所有有效性标记均保留在 环量核查结果 [original local attachment: junction_circulation.json]。

Q>50 正核与 1.5–2.5 等效半径流体环带比较：baseline 的 ΔCp(core−ring) 在 t/T=0 为 +0.0040，t/T=0.5 为 −0.0052。因而不能预设“这个根部核一直是低压核”。此环带比较不是表面压力积分；已有 1.596 的升力比和区域压力收支见原机制报告。

## 受力图的另一处口径问题

用户原图灰线是 baseline−No Fin：L(B+F)−L(B₀)=L(F in baseline)+[L(B in baseline)−L(B₀)]。它含身体载荷因鳍而改变的部分；若还有其他表面，其响应差也在内，不能直接标为独立胸鳍升力。原先报告的 1.596 使用实际胸鳍表面压力载荷，统计口径不同。静态数据的压力升力比约 1.105，也并非严格 1；不要再把静态增益写成精确零。

## 按 Dong 路径进行的对应与新颖性边界

[](https://doi.org/10.1017/jfm.2016.175)

[](https://doi.org/10.1088/1748-3190/ab2208)

[](https://doi.org/10.1017/jfm.2017.533)

[](https://doi.org/10.1103/physrevfluids.8.073101)

[](https://doi.org/10.1063/1.5134689)

[](https://doi.org/10.1155/2016/2721968)

[](https://doi.org/10.3390/biomimetics9110682)

[](https://arxiv.org/html/2409.04147v1)
[](https://arxiv.org/abs/2604.21054)

[](https://arxiv.org/html/2604.21054v1)
[](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/computational-modelling-and-analysis-of-the-hydrodynamics-of-a-highly-deformable-fish-pectoral-fin/CC6284AAF6CD4771F20D4253C54718F3)

| 原始研究 | 获取范围 | 已有描述／分析方法 | 对当前命名的含义 |
| --- | --- | --- | --- |
| Liu, Dong & Li (2016)Vortex dynamics and new lift enhancement mechanism of wing–body interaction in insect forward flight本地原文 [original local attachment: file:///[local path omitted]] | 本地全文；相关正文和图注核对Figs. 6–7, §3.2 | TXV 和 PBV 按身体生成位置及演化命名；WBI 是相互作用机制。WB 与 WN 的根涡本身也发生变化。 | 不能把 Orca 身侧涡直接等同于昆虫 TXV／PBV。 |
| Wang et al. (2019)Computational investigation of wing-body interaction and its lift enhancement effect in hummingbird forward flight本地原文 [original local attachment: file:///[local path omitted]] | 本地全文；相关正文和图注核对Figs. 5–7 | L 形 TXV／nTXV 连接胸部与翼根；用旋向截面、速度和压力核对。WO 中也存在 RV。 | 与当前鳍根桥接结构最接近的 Dong 组对照之一；形似并不足以同名。 |
| Liu et al. (2017)Computational analysis of vortex dynamics and performance enhancement due to body–fin and fin–fin interactions in fish-like locomotion本地原文 [original local attachment: file:///[local path omitted]] | 本地全文；相关正文和图注核对§4.3–4.4, Figs. 13–17 | 跟踪 PBV 与尾鳍 LEV 的相遇，以截面涡量与环量连接压力和载荷。环量边界取峰值涡量的 10% 闭合等值线。 | 应照搬证据链，不应在 Orca 未验证输运来源前照搬 wake capture 结论。 |
| Guo et al. (2023)Vortex dynamics and fin-fin interactions resulting in performance enhancement in fish-like propulsion本地原文 [original local attachment: file:///[local path omitted]] | 本地全文；相关正文和图注核对Abstract; vortex/pressure sections | 用实验与计算压力对照，追踪臀鳍所脱涡与尾鳍相互作用，并分解身体和鳍的收益。 | 不能把任何增强尾迹都当成新涡；需要明确作用对象及受力面。 |
| Li, Dong & Cheng (2020)Tip vortices formation and evolution of rotating wings at low Reynolds numbers | 作者所在大学公开摘要Abstract | 增强 LEV 可以伴随反旋 TEV 和次级尖端涡，形成反旋涡环对。 | 次级涡、反旋对、涡环分支本身不是全新概念。 |
| Li & Su (2016)Fluid Dynamics of Biomimetic Pectoral Fin Propulsion Using Immersed Boundary Method | 公开论文正文Vortex-structure discussion | 胸鳍边缘涡与尖端涡组成连接结构，后续周期有新旧涡环和分支。 | 编号应区分实例／周期；不能从三条弯管推出三类 TEV。 |
| Lin, Gao & Zhang (2024)A Numerical Study of Flow Past a Wall-Mounted Dolphin Dorsal Fin at Low Reynolds Numbers公开全文 | 公开作者预印本全文；期刊元数据§3.1; Figs. 7–9 | 壁面背鳍已有半马蹄、arch、hairpin 及尾迹结构。 | U 形／桥形结构不能仅凭形态主张新颖；该工况为静态壁面背鳍，与 Orca 胸鳍不同。 |
| He & Breuer (2026，预印本)Hydrodynamic loads and vortex evolution from a bio-inspired pectoral fin near a solid body公开全文 | 公开预印本全文；未作为已同行评审论文处理§3.2; Figs. 5–7 | 摆动鳍近身体实验已描述主涡、反旋次级涡、身体反向剪切、角区吸力及射流。 | 这是最直接的新颖性约束之一；近二维端板模型的 tip 不能直接对应三维胸鳍展向尖端。 |
| Dong et al. (2010，基础文献)Computational modelling and analysis of the hydrodynamics of a highly deformable fish pectoral fin | 出版社／作者机构摘要Abstract | 复杂三维胸鳍运动可产生长期附着的强尖端涡；把涡动力学与表面力分布对应。 | 胸鳍涡的命名须考虑实际运动和生成边缘，不能完全套用二维翼型。 |

可采用的中性描述：body-associated reorganization and intensification of the pectoral-fin vortex system。桥接形态、反旋对、次级涡本身均有先例。若后续证实身体运动诱发一个独立分支，可描述为 body-induced secondary junction vortex；但要同时给出来源、相位、旋向和与原涡的连接证据。经典 junction horseshoe vortex 还需要根部上游边界层绕障碍分离的证据；U 形外观不满足这一条件。vortex reconnection 则需要涡线／拓扑演化证据，Q 壳相接不够。

现有数据的未闭合处是：一个周期、Δt/T=0.0625 的时间分辨率；原网格不同且未做局部涡核网格收敛；两层壁面排除区遮挡涡量源；暂无跨相位材料追踪与全涡管端点辨认。因此报告既不宣称发现新涡型，也不宣称已经证明没有额外生成的涡。相比盲目增加 Q 图片，最有判别力的后续数据是在鳍根涡形成阶段保存更密的三维速度／压力及壁面附近流场，并按有／无身体运动的控制对照定位生成源。

## 检索覆盖与复核清单

本地源：research literacture/04_目录与记录/authors/Haibo_Dong/全部论文清单.csv，经 manifest 映射现存 PDF。2016–2026 年 196 条目录记录，49 篇有可提取全文；对这些全文做逐页关键词筛查，再精读直接相关研究。包含会议记录和不同发表阶段，不能按 196 篇独立期刊论文计。2016 年较早的关键翼身论文作为约十年窗口的补充保留；另有 2010 基础文献。

2026-10-22 的公开检索组合涵盖 wing-body interaction、root/thorax/posterior-body vortex、fin-body junction、pectoral fin secondary vortex、wall-mounted dolphin fin、tip vortex loops、LEV/TEV origin 等；只用原始论文、作者预印本、出版社和作者机构页面支撑判断。2020 Dong 尖端涡论文补得摘要；2026 He–Breuer 按预印本处理。对缺失的 147 条本地全文没有假装已阅读全文；检索也不是对全领域的穷尽性新颖性证明。

展开全部 196 条目录记录及本地全文状态

| 年 | 题名 | 状态 | 目录 ID |
| --- | --- | --- | --- |
| 2026 | Design of A Lightweight Robotic Tensegrity Morphing Airfoil | 仅目录记录；缺本地全文 | W7125921905 |
| 2026 | Fins in formation: hydrodynamic impact of median fins in in-line fish swimming | 全文提取／术语筛查 | W7119758505 |
| 2026 | Hydrodynamic interactions of low-aspect-ratio oscillating panels in a tip-to-tip formation | 全文提取／术语筛查 | W7133213453 |
| 2026 | Minimal-actuation feather star–inspired soft swimmers for multimodal 3D maneuverability | 仅目录记录；缺本地全文 | W7220864324 |
| 2026 | Passive thrust production in animal propulsion: exploitation of velocity gradients for energy conservation | 仅目录记录；缺本地全文 | W7213412583 |
| 2026 | Tuna-like swimmers experience a fluid-mediated stable side-by-side formation | 仅目录记录；缺本地全文 | W7166545250 |
| 2025 | A Computational Study of In-Phase and Anti-Phase Interactions of Fish in a Phalanx School | 仅目录记录；缺本地全文 | W4408079097 |
| 2025 | An Untethered Bioinspired Robotic Tensegrity Dolphin with Multi-Flexibility Design for Aquatic Locomotion | 全文提取／术语筛查 | W4411019178 |
| 2025 | Computational analysis of fish-foil pairing and wake energy extraction in low-speed flow | 全文提取／术语筛查 | W4414114672 |
| 2025 | Fish schools in a vertical diamond formation: Effect of vertical spacing on hydrodynamic interactions | 全文提取／术语筛查 | W4409766453 |
| 2025 | Hydrodynamic Interactions Between Schooling Fish in a Three-Dimensional Non-Planar Cluster | 全文提取／术语筛查 | W4414461969 |
| 2025 | Numerical Analysis of the Aerodynamic Interactions in Tandem Flying Snake Airfoils | 仅目录记录；缺本地全文 | W4408343735 |
| 2025 | Physics-Informed Neural Networks Simulation and Validation of Airflows in Three-Dimensional Upper Respiratory Tracts | 全文提取／术语筛查 | W4416618243 |
| 2025 | Tuna-Like Swimmers Experience a Fluid-Mediated Stable Side-by-Side Formation | 全文提取／术语筛查 | W7116944321 |
| 2025 | Tuna-Like Swimmers Experience a Fluid-Mediated Stable Side-by-Side Formation | 仅目录记录；缺本地全文 | W7117078958 |
| 2025 | Vortex dynamics in wake–body and wake-fin interactions of tuna-like staggered swimming | 全文提取／术语筛查 | W4411986689 |
| 2025 | Withdrawn: Design of a Lightweight Robotic Tensegrity Morphing Airfoil | 仅目录记录；缺本地全文 | W4407415664 |
| 2024 | Computational Analysis of Foil-Fish Interactions in Low Reynolds Number Flow | 仅目录记录；缺本地全文 | W4402412913 |
| 2024 | Computational Study of In-Phase and Anti-Phase Interactions of Fish in a Phalanx School | 全文提取／术语筛查 | W4402413939 |
| 2024 | Effects of body shape on hydrodynamic interactions in a dense diamond fish school | 全文提取／术语筛查 | W4392715753 |
| 2024 | Flow Analysis of Central Venous Outflow Tract: A New Approach to Understanding Pulse‐Synchronous Tinnitus | 仅目录记录；缺本地全文 | W4400495767 |
| 2024 | Hydrodynamic Interactions in Fish-Like Robotic Swarms With Flexible Propulsors | 全文提取／术语筛查 | W4402413881 |
| 2024 | Hydrodynamic analysis of fin–fin interactions in two-manta-ray schooling in the vertical plane | 仅目录记录；缺本地全文 | W4390574153 |
| 2024 | Hydrodynamics of Active Fin Control in Fish-Like Swimming | 全文提取／术语筛查 | W4406763179 |
| 2024 | Morphokinematic Impacts on Hydrodynamics and Performance of Tuna-Inspired Robots | 仅目录记录；缺本地全文 | W4406763030 |
| 2024 | Patient-specific computational fluid dynamic modeling accurately predicts systolic pressure gradients across stenotic right ventricular outflow tracts | 仅目录记录；缺本地全文 | W4401475420 |
| 2024 | Postural Change of the Annual Cicada (Tibicen linnei) Helps Facilitate Backward Flight | 全文提取／术语筛查 | W4394821619 |
| 2024 | Spontaneous snapping-induced jet flows for fast, maneuverable surface and underwater soft flapping swimmer | 全文提取／术语筛查 | W4405023162 |
| 2024 | Unraveling hydrodynamic interactions in fish schools: A three-dimensional computational study of in-line and side-by-side configurations | 全文提取／术语筛查 | W4401718563 |
| 2024 | Vertical bending and aerodynamic performance in flying snake-inspired aerial undulation | 仅目录记录；缺本地全文 | W4404311433 |
| 2024 | Wing deformation improves aerodynamic performance of forward flight of bluebottle flies flying in a flight mill | 仅目录记录；缺本地全文 | W4400738091 |
| 2024 | Wing–antenna interaction reduces odour fatigue in butterfly odour-tracking flight | 全文提取／术语筛查 | W4404020881 |
| 2023 | A Multiphysics Approach to Understanding Chemoreception in Bio-Robotic Fish Schools | 仅目录记录；缺本地全文 | W4391551153 |
| 2023 | Aerodynamic Forces and Wake Analysis of Wing Damaged Flapping Flight | 仅目录记录；缺本地全文 | W4317597214 |
| 2023 | Computation study about the interaction between the tandem flying snake airfoils with dynamic motion | 仅目录记录；缺本地全文 | W4317651259 |
| 2023 | Convergence in Gliding Animals: Morphology, Behavior, and Mechanics | 仅目录记录；缺本地全文 | W4322749678 |
| 2023 | Data-Driven Discovery of Anomaly-Sensitive Parameters from Uvula Wake Flows Using Wavelet Analyses and Poincaré Maps | 全文提取／术语筛查 | W4388226621 |
| 2023 | Direct numerical simulations and flow-pressure acoustic analyses of flapping-uvula-induced flow evolutions within normal and constricted pharynx | 全文提取／术语筛查 | W4324144837 |
| 2023 | Energy harvesting of passively pitching foil through vortex capturing | 仅目录记录；缺本地全文 | W4317611154 |
| 2023 | Flow Interactions and Wake Structures in Large Groups of Undulating Foils | 仅目录记录；缺本地全文 | W4317632544 |
| 2023 | Geometric Characteristics of Flapping Foils for Enhanced Propulsive Efficiency | 仅目录记录；缺本地全文 | W4322492662 |
| 2023 | Hydrodynamic Performance and Wake Topology of Schooling Fish in Three-Dimensional Formations | 仅目录记录；缺本地全文 | W4323865006 |
| 2023 | Hydrodynamics of Metachronal Rowing at Intermediate Reynolds Numbers | 仅目录记录；缺本地全文 | W4391551053 |
| 2023 | The Hydrodynamic Benefits Achieved Through Lateral Spacing of Schooling Manta Rays | 仅目录记录；缺本地全文 | W4391551001 |
| 2023 | Thrust generation and propulsive efficiency in dolphin-like swimming propulsion | 全文提取／术语筛查 | W4383312750 |
| 2023 | Vortex dynamics and fin-fin interactions resulting in performance enhancement in fish-like propulsion | 全文提取／术语筛查 | W4385067820 |
| 2022 | Aerodynamics and Flow Characteristics of Flying Snake Gliding with Undulating Locomotion | 仅目录记录；缺本地全文 | W4206576942 |
| 2022 | An Efficient Tree-Topological Local Mesh Refinement on Cartesian Grids for Multiple Moving Objects in Incompressible Flow | 仅目录记录；缺本地全文 | W4319075154 |
| 2022 | Bio-Inspired Propulsion: Towards Understanding the Role of Pectoral Fin Kinematics in Manta-like Swimming | 全文提取／术语筛查 | W4223999777 |
| 2022 | Body Shape Effects on the Hydrodynamic Performance of Bio-Inspired Undulating Swimmers | 仅目录记录；缺本地全文 | W4312506599 |
| 2022 | Computational Modeling and Hydrodynamic Analysis of Fish Schools in Three-Dimensional Arrangements | 仅目录记录；缺本地全文 | W4312477586 |
| 2022 | Computational analysis of vortex dynamics and aerodynamic performance in flying-snake-like gliding flight with horizontal undulation | 仅目录记录；缺本地全文 | W4308145921 |
| 2022 | Computational study on the gliding flight of a damselfly | 仅目录记录；缺本地全文 | W4206193907 |
| 2022 | Effect of Tunable Stiffness on the Hydrodynamic Performance of a Tuna Tail Informed Flexible Propulsor | 全文提取／术语筛查 | W4205749095 |
| 2022 | Effects of phase difference on hydrodynamic interactions and wake patterns in high-density fish schools | 全文提取／术语筛查 | W4307983670 |
| 2022 | Flow Dynamics and Acoustics from Glottal Vibrations at Different Frequencies | 全文提取／术语筛查 | W4307944114 |
| 2022 | Links between language, fluid dynamics, and the airflows that transport pathogens | 仅目录记录；缺本地全文 | W4312497899 |
| 2022 | Modeling and Computation of Batoid Swimming Inspired Pitching Impact on Wake Structure and Hydrodynamic Performance | 全文提取／术语筛查 | W4313158343 |
| 2022 | Role of the caudal peduncle in a fish-inspired robotic model: how changing stiffness and angle of attack affects swimming performance | 仅目录记录；缺本地全文 | W4303474671 |
| 2022 | Sound generated by flow over two traveling wavy foils in a side-by-side arrangement | 仅目录记录；缺本地全文 | W4310168494 |
| 2022 | Wing Kinematics and Unsteady Aerodynamics of a Hummingbird Pure Yawing Maneuver | 全文提取／术语筛查 | W4292547544 |
| 2021 | A Versatile IBM-Based AMR Method for Studying Human Snoring | 仅目录记录；缺本地全文 | W3202049801 |
| 2021 | A sharp interface immersed boundary method for flow-induced noise prediction using acoustic perturbation equations | 仅目录记录；缺本地全文 | W3171495233 |
| 2021 | A study of three-dimensional maneuvering hydrodynamics inspired by a larval zebrafish model | 仅目录记录；缺本地全文 | W3186369830 |
| 2021 | An Improved Level-Set-Based Immersed Boundary Reconstruction Method for Computing Bio-Inspired Underwater Propulsion | 仅目录记录；缺本地全文 | W3202423071 |
| 2021 | Computational Investigation of Thrust Production of a Dolphin at Various Swimming Speeds | 仅目录记录；缺本地全文 | W3203581578 |
| 2021 | Correction: A study of three-dimensional maneuvering hydrodynamics inspired by a larval zebrafish model | 仅目录记录；缺本地全文 | W3193436803 |
| 2021 | Effect of Formation Pattern on Schooling Energetics in Fish-Like Swimming | 仅目录记录；缺本地全文 | W4251198940 |
| 2021 | Effects of Chordwise or Spanwise Non-uniform Stiffness on Propulsive Performance of Square Foil | 仅目录记录；缺本地全文 | W3119522939 |
| 2021 | Larger wavelengths suit hydrodynamics of carangiform swimmers | 全文提取／术语筛查 | W3092451939 |
| 2021 | Leveraging statistical shape modeling in computational respiratory dynamics: Nanomedicine delivery in remodeled airways | 仅目录记录；缺本地全文 | W3149146005 |
| 2021 | Propulsive performance and vortex wakes of multiple tandem foils pitching in-line | 全文提取／术语筛查 | W3214154487 |
| 2021 | Video: Effects of wing-induced flow on the odor plume structures in an upwind surging flight of monarch butterfly | 仅目录记录；缺本地全文 | W3205372381 |
| 2021 | Video: Swimming in Schools: Computational Analysis of Hydrodynamic Interactions in Trout Schools | 仅目录记录；缺本地全文 | W3215234050 |
| 2021 | Wake Structures and Effect of Hydrofoil Shapes in Efficient Flapping Propulsion | 仅目录记录；缺本地全文 | W3204872566 |
| 2021 | Why do anguilliform swimmers perform undulation with wavelengths shorter than their bodylengths? | 仅目录记录；缺本地全文 | W3139361763 |
| 2020 | A CFD simulation of a flying snake gliding with body undulations | 仅目录记录；缺本地全文 | W3119905086 |
| 2020 | A Versatile Level Set Based Immersed Boundary Reconstruction for Bio-Inspired Flow Applications | 仅目录记录；缺本地全文 | W2996983856 |
| 2020 | A parallel Local-Adaptive-Mesh-Refinement-enabled Immersed Boundary Method for biological flows | 仅目录记录；缺本地全文 | W3108400629 |
| 2020 | Aerodynamics and Wing-Wing Interaction during the Pre-Ovipository Flight of the Damselfly | 仅目录记录；缺本地全文 | W2998328613 |
| 2020 | Airflow Characterization Combined With 3D Reconstruction of Rabbit Vocal Cord Based on In-Vivo Phonatory Experiments | 仅目录记录；缺本地全文 | W3107534260 |
| 2020 | Alveolar size effects on nanoparticle deposition in rhythmically expanding-contracting terminal alveolar models | 仅目录记录；缺本地全文 | W3021510721 |
| 2020 | CFD Analysis of Micro-Channel Flows for Applications in Environmental Structure Cracking | 仅目录记录；缺本地全文 | W2997513337 |
| 2020 | Computational analysis of hydrodynamic interactions in a high-density fish school | 全文提取／术语筛查 | W3111867573 |
| 2020 | Dual functions of insect wings in an odor-guided aeronautic navigation | 仅目录记录；缺本地全文 | W3004647005 |
| 2020 | Effect of Formation Pattern on Schooling Energetics in Fish-Like Swimming | 仅目录记录；缺本地全文 | W3118317606 |
| 2020 | Effect of Head Posture on Airflow and Pressure Behavior of Human Upper Airway During Snoring | 仅目录记录；缺本地全文 | W3180882564 |
| 2020 | Effects of Chordwise Non-uniform Stiffness on Propulsive Performance of Passively-flexing Square Foil | 仅目录记录；缺本地全文 | W2997284667 |
| 2020 | Effects of Intermittent Swimming Gait in Fish-Like Locomotion | 仅目录记录；缺本地全文 | W2998233015 |
| 2020 | Effects of Torsional Stiffness on Hydrodynamic Performance of a Heaving and Passive Pitching Trapezoidal Panel | 仅目录记录；缺本地全文 | W2997443052 |
| 2020 | Effects of tunable stiffness on the hydrodynamics and flow features of a passive pitching panel | 仅目录记录；缺本地全文 | W3094644696 |
| 2020 | Extracting signature responses from respiratory flows: Low‐dimensional analyses on Direct Numerical Simulation‐predicted wakes of a flapping uvula | 仅目录记录；缺本地全文 | W3092884284 |
| 2020 | Flowtaxis in the wakes of oscillating airfoils | 全文提取／术语筛查 | W3048004943 |
| 2020 | Hydrodynamics of median-fin interactions in fish-like locomotion: Effects of fin shape and movement | 全文提取／术语筛查 | W2999224801 |
| 2020 | Kinematics and Hydrodynamics of a Dolphin in Forward Swimming | 仅目录记录；缺本地全文 | W3035431155 |
| 2020 | Modal decompositions of the kinematics of Crevalle jack and the fluid–caudal fin interaction | 全文提取／术语筛查 | W3093221717 |
| 2020 | Numerical simulation of acoustic scattering and flow-induced noise with sharp interface immersed boundary method | 仅目录记录；缺本地全文 | W2998587199 |
| 2020 | Numerical study of hydrodynamics interactions in a diamond fish school | 仅目录记录；缺本地全文 | W3120829194 |
| 2020 | Scaling the vorticity dynamics in the leading-edge vortices of revolving wings with two directional length scales | 仅目录记录；缺本地全文 | W3117323563 |
| 2020 | The reverse flight of a monarch butterfly (Danaus plexippus) is characterized by a weight-supporting upstroke and postural changes | 仅目录记录；缺本地全文 | W3037930861 |
| 2020 | Tip vortices formation and evolution of rotating wings at low Reynolds numbers | 仅目录记录；缺本地全文 | W3006217317 |
| 2020 | Tuna Swimming: An Enhanced Hydrodynamic Performance by Finlets | 仅目录记录；缺本地全文 | W3112818703 |
| 2020 | Tuna locomotion: a computational hydrodynamic analysis of finlet function | 全文提取／术语筛查 | W3014515713 |
| 2020 | Twist Morphing Effect on Propulsive Performance of Bio-Inspired Pitching-Rolling Plates | 全文提取／术语筛查 | W2996844467 |
| 2019 | An image-guided computational approach to inversely determine in vivo material properties and model flow-structure interactions of fish fins | 仅目录记录；缺本地全文 | W2942527978 |
| 2019 | Comparison of Geometric Parameterization Methods for Optimal Shape Design in Efficient Flapping Propulsion | 仅目录记录；缺本地全文 | W2949148518 |
| 2019 | Computational analysis of a flapping uvula on aerodynamics and pharyngeal wall collapsibility in sleep apnea | 仅目录记录；缺本地全文 | W2963072905 |
| 2019 | Computational investigation of wing-body interaction and its lift enhancement effect in hummingbird forward flight | 全文提取／术语筛查 | W2944837325 |
| 2019 | Computational study of fish-shaped panel with simultaneously heaving and bending motion | 仅目录记录；缺本地全文 | W2909535718 |
| 2019 | Effect of Diameter Size on the Unsteady Characteristics of Pulsatile Flow in Curved Human Blood Vessel | 仅目录记录；缺本地全文 | W2991008459 |
| 2019 | Effect of Joint Stiffness on Propulsive Performance of Fish-Shaped Panel | 仅目录记录；缺本地全文 | W2950396708 |
| 2019 | Effect of Uvula Length on Airflow and Pressure Oscillation in a Human Pharynx Model | 仅目录记录；缺本地全文 | W3003003786 |
| 2019 | Foil shapes for efficient fish-like propulsion | 仅目录记录；缺本地全文 | W2910292385 |
| 2019 | How dorsal fin sharpness affects swimming speed and economy | 全文提取／术语筛查 | W2972951149 |
| 2019 | Hydrodynamics and Flow Characterization of Tuna-Inspired Propulsion in Forward Swimming | 仅目录记录；缺本地全文 | W2990762720 |
| 2019 | Introduction to the Special Issue on Advances in Bio-Inspired Propulsion | 仅目录记录；缺本地全文 | W2971642262 |
| 2019 | Large-Amplitude Intermittent Swimming in Viscous and Inviscid Flows | 仅目录记录；缺本地全文 | W2948941789 |
| 2019 | Numerical Investigation on Hydrodynamic Performance of Flapping Plates with Non-Uniform Spanwise Flexibility Using Fluid Structure Interaction | 仅目录记录；缺本地全文 | W2953385010 |
| 2019 | Passive Pitching Mechanism of Three-Dimensional Flapping Wings in Hovering Flight | 仅目录记录；缺本地全文 | W2991101844 |
| 2019 | Radial planetary vorticity tilting in the leading-edge vortex of revolving wings | 仅目录记录；缺本地全文 | W2941692321 |
| 2019 | Study on the passive pitching mechanism of different forms of flapping motion in turning flight | 仅目录记录；缺本地全文 | W2950431114 |
| 2019 | Ventilation Modulation and Nanoparticle Deposition in Respiratory and Olfactory Regions of Rabbit Nose | 全文提取／术语筛查 | W2992778505 |
| 2019 | Wake Structure and Aerodynamic performance of Passively Pitching Revolving Plates | 仅目录记录；缺本地全文 | W2908863984 |
| 2018 | A Numerical Study of the Sound and Force Production of Flexible Insect Wings | 全文提取／术语筛查 | W2898690810 |
| 2018 | A balance between aerodynamic and olfactory performance during flight in Drosophila | 仅目录记录；缺本地全文 | W2885122802 |
| 2018 | A level set based boundary reconstruction method for 3-D bio-inspired flow simulations with sharp-interface immersed boundary method | 仅目录记录；缺本地全文 | W4232977709 |
| 2018 | Adjoint-Based Optimization for Thrust Performance of Three-Dimensional Pitching–Rolling Plate | 仅目录记录；缺本地全文 | W2892255035 |
| 2018 | Adjoint-based optimization for thrust performance of a three-dimensional pitching-rolling plate | 全文提取／术语筛查 | W2950781753 |
| 2018 | An Integrated High-fidelity Approach for Modeling Flow-structure Interaction in Biological Propulsion and its Strong Validation | 仅目录记录；缺本地全文 | W2784054033 |
| 2018 | Bifurcations and route to chaos for flow over an oscillating airfoil | 仅目录记录；缺本地全文 | W2802996015 |
| 2018 | Computational Analysis on Aerodynamics and Vortex Formation of Sleep Apnea | 仅目录记录；缺本地全文 | W2898554297 |
| 2018 | Computational Study of the Three-Dimensional Wake and Performance of a Trapezoidal Pitching Panel | 仅目录记录；缺本地全文 | W2783731496 |
| 2018 | Correction: A level set based boundary reconstruction method for 3-D bio-inspired flow simulations with sharp-interface immersed boundary method | 仅目录记录；缺本地全文 | W2810284067 |
| 2018 | Correction: Effects of flapping waveforms on the performance of intermittent swimming in viscous flows | 仅目录记录；缺本地全文 | W2791764055 |
| 2018 | Dual functions of insect wings: balancing aerodynamics and olfaction | 仅目录记录；缺本地全文 | W2922364795 |
| 2018 | Effects of flapping waveforms on the performance of intermittent swimming in viscous flows | 全文提取／术语筛查 | W2783810442 |
| 2018 | Effects of glottis motion on airflow and energy expenditure in a human upper airway model | 仅目录记录；缺本地全文 | W2800675436 |
| 2018 | Effects of uvula morphology on the fluid dynamics in human snoring | 仅目录记录；缺本地全文 | W3202364971 |
| 2018 | Fast level-set-based boundary reconstruction of highly complex body configurations for studying bio-inspired locomotion using a sharp-interface immersed boundary method | 仅目录记录；缺本地全文 | W2911929617 |
| 2018 | Flying in reverse: kinematics and aerodynamics of a dragonfly in backward free flight | 仅目录记录；缺本地全文 | W2809787148 |
| 2018 | How dorsal fin sharpness affects swimming speed and efficiency | 仅目录记录；缺本地全文 | W2922824445 |
| 2018 | Kinematics and Hydrodynamics of Invertebrate Undulatory Swimming | 仅目录记录；缺本地全文 | W2898562827 |
| 2018 | On the high-performance swimming of a tuna-inspired underwater vehicle Part II: Computational analysis | 仅目录记录；缺本地全文 | W2918054485 |
| 2018 | On the hydrodynamics and nonlinear interaction between fish in tandem configuration | 仅目录记录；缺本地全文 | W2794939388 |
| 2018 | Optimizing foil shape for fish-like propulsion Part I: Experimental performance | 仅目录记录；缺本地全文 | W2911341507 |
| 2018 | Optimizing foil shape for fish-like propulsion Part II: Computational optimization and hydrodynamics | 仅目录记录；缺本地全文 | W2916309783 |
| 2018 | Passive Pitching of a Flapping Wing in Turning Flight | 仅目录记录；缺本地全文 | W2793969331 |
| 2018 | Phonatory experiments and computational analysis of human larynges using ex vivo models | 仅目录记录；缺本地全文 | W2908765469 |
| 2018 | Video: Gliding with body undulations: 3D flow simulations of a flying snake | 仅目录记录；缺本地全文 | W2908680189 |
| 2018 | Vorticity dynamics of revolving wings: The role of planetary vortex tilting on the stability of the leading-edge vortex | 全文提取／术语筛查 | W2774615069 |
| 2018 | Wake capture and aerodynamics of passively pitching tandem flapping plates | 仅目录记录；缺本地全文 | W2810757242 |
| 2017 | A highly efficient sharp-interface immersed boundary method with adaptive mesh refinement for bio-inspired flow simulations. | 仅目录记录；缺本地全文 | W2775217251 |
| 2017 | Aerodynamics and flow features of a damselfly in takeoff flight | 仅目录记录；缺本地全文 | W2734565635 |
| 2017 | An Adjoint-Based Approach to Study a Flexible Flapping Wing in Pitching-Rolling Motion | 仅目录记录；缺本地全文 | W2774182650 |
| 2017 | Computational analysis of vortex dynamics and performance enhancement due to body–fin and fin–fin interactions in fish-like locomotion | 全文提取／术语筛查 | W2754658223 |
| 2017 | Computational modeling and validation of human nasal airflow under various breathing conditions | 仅目录记录；缺本地全文 | W2753810598 |
| 2017 | Effects of Surface Morphing on the Wake Structure and Performance of Flapping Plates | 仅目录记录；缺本地全文 | W2620960784 |
| 2017 | Effects of aspect ratio and angle of attack on tip vortex structures and aerodynamic performance for rotating flat plates | 仅目录记录；缺本地全文 | W2621264214 |
| 2017 | Enhancing the Efficiency of Bio-Inspired Propulsion via Intermittent Swimming Gaits | 全文提取／术语筛查 | W2620594989 |
| 2017 | Flow-structure Interaction Modeling of a Fish Caudal Fin during Steady Swimming | 仅目录记录；缺本地全文 | W2771650710 |
| 2017 | How Insects Initiate Flight: Computational Analysis of a Damselfly in Takeoff Flight | 仅目录记录；缺本地全文 | W2774544932 |
| 2017 | How can a source be located by responding to local information in its hydrodynamic trail? | 仅目录记录；缺本地全文 | W2782931905 |
| 2017 | How swimming near a curved body could improve bio-inspired propulsion | 仅目录记录；缺本地全文 | W2771141830 |
| 2017 | Kinematics and Hydrodynamics of Burst-and-Coast Strategy in Carangiform Swimming | 仅目录记录；缺本地全文 | W2774766008 |
| 2017 | Kinematics and Hydrodynamics of Mobuliform Swimming: Oscillatory Winged Propulsion by Large Pelagic Batoids | 仅目录记录；缺本地全文 | W2767701619 |
| 2017 | Low Dimensional kinematic models and hydrodynamic performance of a Trout in steady swimming | 仅目录记录；缺本地全文 | W2775159854 |
| 2017 | Nonlinear Hydrodynamic Interaction in Fish Swimming | 仅目录记录；缺本地全文 | W2568151833 |
| 2017 | Numerical study of dynamic glottis and tidal breathing on respiratory sounds in a human upper airway model | 仅目录记录；缺本地全文 | W2766479548 |
| 2017 | On the passive pitching mechanism in turning flapping flights using a torsional spring model | 仅目录记录；缺本地全文 | W2620875831 |
| 2017 | On the sound production and radiation of bio-inspired propulsors | 仅目录记录；缺本地全文 | W2621341207 |
| 2017 | Optimized Body Deformation in Dragonfly Maneuvers | 全文提取／术语筛查 | W2737243004 |
| 2017 | Proper Orthogonal Decomposition Analysis of Flapping Hovering Wings | 仅目录记录；缺本地全文 | W2569197441 |
| 2017 | Proper Orthogonal Decomposition Analysis on the Hydrodynamic Effects of the Propulsor Geometries in Flapping Propulsion | 全文提取／术语筛查 | W2569039375 |
| 2017 | Propulsive Performance of Pitching Panels with Bio-Inspired Passive Directional Flexibility | 仅目录记录；缺本地全文 | W2620828093 |
| 2017 | The effect of wing flexibility on sound generation of flapping wings | 仅目录记录；缺本地全文 | W2744164699 |
| 2017 | Wing kinematics measurement and aerodynamics of a dragonfly in turning flight | 仅目录记录；缺本地全文 | W2569694275 |
| 2016 | Cetacean Swimming with Prosthetic Limbs | 仅目录记录；缺本地全文 | W2568099681 |
| 2016 | Cicada (Tibicen linnei) steers by force vectoring | 仅目录记录；缺本地全文 | W2256128762 |
| 2016 | Computational Analysis of 3D Fin-Fin Interaction in Fish’s Steady Swimming | 仅目录记录；缺本地全文 | W2560832327 |
| 2016 | Effect of body configuration and deformation on the dynamic stability of flapping flight | 仅目录记录；缺本地全文 | W2408777171 |
| 2016 | Effects of Tail Geometries on the Performance and Wake Pattern in Flapping Propulsion | 仅目录记录；缺本地全文 | W2560011625 |
| 2016 | Fin-Body Interaction and its Hydrodynamic Benefits in Fish's Steady Swimming | 仅目录记录；缺本地全文 | W2563637393 |
| 2016 | Hydrodynamic Performance of Aquatic Flapping: Efficiency of Underwater Flight in the Manta | 全文提取／术语筛查 | W2473194771 |
| 2016 | Hydrodynamics of a tandem fish school with asynchronous undulation of individuals | 仅目录记录；缺本地全文 | W2483612304 |
| 2016 | Investigation of Body-involved Lift Enhancement in Bio-inspired Flapping Flight | 仅目录记录；缺本地全文 | W2564833347 |
| 2016 | Low-dimensional Modeling and Aerodynamics of Flexible Wings in Flapping Flight | 仅目录记录；缺本地全文 | W2409560335 |
| 2016 | New Insights on Insect's Silent Flight. Part I: Vortex Dynamics and Wing Morphing | 仅目录记录；缺本地全文 | W2560911129 |
| 2016 | New insights into insect's silent flight. Part II: sound source and noise control | 仅目录记录；缺本地全文 | W2609199184 |
| 2016 | Physics-based Morphology Analysis and Adjoint Optimization of Flexible Flapping Wings | 仅目录记录；缺本地全文 | W2763192366 |
| 2016 | Proper Orthogonal Decomposition Analysis of 3-D Wake Structures in a Pitching-Rolling Plate | 仅目录记录；缺本地全文 | W2325334706 |
| 2016 | Quantification and Analysis of Propulsive Wake Topologies in Finite Aspect-Ratio Pitching-Rolling Plates | 仅目录记录；缺本地全文 | W2415533830 |
| 2016 | Quantification of wing and body kinematics in connection to torque generation during damselfly yaw turn | 仅目录记录；缺本地全文 | W2553871379 |
| 2016 | Turning on a dime: Asymmetric vortex formation in hummingbird maneuvering flight | 全文提取／术语筛查 | W2519423505 |
| 2016 | Video: Swimming Upstream: Computational Hydrodynamics of Trout Locomotion | 仅目录记录；缺本地全文 | W2554154107 |
| 2016 | Vortex dynamics and new lift enhancement mechanism of wing–body interaction in insect forward flight | 全文提取／术语筛查 | W2338466152 |
| 2016 | ベクトル力による去勢牛セミ(Tibicen linnei)【Powered by NICT】 | 仅目录记录；缺本地全文 | W2744782818 |

完整筛查 CSV [original local attachment: literature/dong_2016_2026_screening.csv] · 关键文献与定位信息 [original local attachment: literature/curated_evidence.json] · 逐相位区域数据 [original local attachment: junction_region.json]

辅助诊断图（非主结果）

以下用于复核坐标、旋向和匹配，不额外支撑“新涡型”的主张。

[Figure unavailable in source: figures/geometric_phase_check.png]

[Figure unavailable in source: figures/section_core_tracks.png]

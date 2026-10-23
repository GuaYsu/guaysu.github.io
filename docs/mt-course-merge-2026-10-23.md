# 大地电磁教程合并记录（2026-10-23）

## 目标与范围

- 合并前：`_posts` 中 60 篇课程相关 Markdown，包括主目录、7 篇专题导读、47 篇逐步教程及 5 篇早期留存反演教程。
- 合并后：保留 1 篇课程主目录，形成 **17 篇正式合并讲义**。
- 完整保留每篇原教程的正文、图片引用、数学公式、例题、练习与概念说明。合并只调整原有 Markdown 的标题层级、跨原章节导航链接与重复的 MathJax include。
- 共 59 篇被替换的源文件以**未修改的 Git Blob**放在 `_mt-course-archive/`，并保留版本历史。
- 59 个旧 permalink 通过 `mt-course-redirects/` 下的普通页面与 `jekyll-redirect-from` 插件跳转到对应合并篇的锚点。重定向页面放在普通页面目录，而不是 `_posts` 中，避免继续占据博客文章列表。

## 合并后 17 篇正式讲义

| 顺序 | 新文件后缀 | 来源章节后缀 |
|---|---|---|
| 1 | `01-foundations` | `01-field-wave`, `01b-trigonometry`, `01c-complex` |
| 2 | `02a-maxwell-math` | `02-maxwell-forward`, `02a-derivatives`, `02b-curl`, `02c-1d-ode` |
| 3 | `02b-layered-earth` | `02d-impedance`, `02e-layers` |
| 4 | `03a-numerical-grid` | `03-finite-difference`, `03a-grid`, `03b-differences`, `03c-nonuniform`, `03d-matrix` |
| 5 | `03b-numerical-forward` | `03e-boundaries`, `03f-multidimensional`, `03g-sensitivity` |
| 6 | `04-field-observation` | `04-field`, `04a-measure`, `04b-station` |
| 7 | `05a-fourier-sampling` | `05-processing`, `05a-fourier`, `05b-sampling`, `05c-windowing` |
| 8 | `05b-impedance-estimation` | `05d-least-squares`, `05e-robust` |
| 9 | `06a-tensor-dimensions` | `06-analysis`, `06a-tensor`, `06b-dimensionality`, `06c-rotation` |
| 10 | `06b-tensor-interpretation` | `06d-invariants`, `06e-distortion`, `06f-phase-tensor` |
| 11 | `06c-strike-static` | `06g-strike`, `06h-static-shift` |
| 12 | `07a-inverse-theory` | `07-inversion`, `07a-inverse-problem`, `07b-nullspace`, `07c-misfit` |
| 13 | `07b-optimization` | `07d-jacobian`, `07e-newton`, `07f-gauss-newton` |
| 14 | `07c-occam-resolution` | `07g-occam`, `07h-validation` |
| 15 | `08a-rock-physics` | `08-interpretation`, `08a-conduction`, `08b-archie`, `08c-global-arrays` |
| 16 | `08b-continental-cases` | `08d-cascadia`, `08e-northamerica`, `08f-sinoprobe` |
| 17 | `08c-field-cases` | `08g-junggar-amt`, `08h-marine`, `08i-evidence` |

源文件 `07a-inverse`, `07b-misfit`, `07c-regularization`, `07e-gauss-newton`, `07f-diagnostics` 是第一版反演教程，未混入当前正式顺序；它们保存在原文归档里，旧链接分别指向含相同主题的扩充版章节。

## 静态检查记录

已检查：
1. Git 树中课程 `_posts` 为 **18 篇**（17 个正文合并篇 + 1 个总目录），不存在重复旧课程短文作为博客文章；
2. `_mt-course-archive/` 中有 59 个原始源文件和 README；
3. `mt-course-redirects/` 中有 59 个入口；
4. 17 篇合并文章都具有原章节 `sec-...` 锚点；
5. 每篇只保留一个 MathJax 加载入口，未发现奇数个 `$$` 展示公式定界符；
6. Jekyll 配置中的 `jekyll-redirect-from` 插件已启用。

待运行时确认：GitHub Pages 构建成功后抽查旧 URL 的实际 HTTP 跳转与移动端文章内锚点滚动。静态检查不能替代在线部署验证。

[查看当前课程总目录](/posts/2026/10/mt-course/)。

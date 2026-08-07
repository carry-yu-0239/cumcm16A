# 模型说明索引

在本文件中汇总各子问题已经确认的模型说明，并链接至相应的模型卡、代码、输出和论文位置。

## 子问题一：静水风载下的刚体—悬链线平衡

- 模型卡：[model_cards/q1.md](../model_cards/q1.md)
- MATLAB 实现：[src/q1/q1_static_model.m](../src/q1/q1_static_model.m)
- 最终表格：[outputs/q1/tables/q1_static_results.csv](../outputs/q1/tables/q1_static_results.csv)
- 曲线与链形：[outputs/q1/figures/](../outputs/q1/figures/)
- 核验日志：[outputs/q1/logs/q1_static_checks.txt](../outputs/q1/logs/q1_static_checks.txt)

以浮标吃水 \(\delta\) 为自变量，浮标风力为
\[
F_w(\delta)=0.625\,[2(2-\delta)]v^2.
\]
对任一均质密闭圆柱构件，取其上端张力的水平分量为 \(H=F_w\)、竖向分量为 \(V_{\rm top}\)，水中有效重力为 \(W\)，由端点力平衡和绕上端的力矩平衡得
\[
V_{\rm bottom}=V_{\rm top}-W,\qquad
\tan\varphi=\frac{H}{V_{\rm top}-W/2}.
\]
钢桶—锚链节点处的重物球作为竖直向下分支载荷；II 型锚链有效单位重力采用
\[
q=\mu g\left(1-\frac{\rho_{\rm water}}{\rho_{\rm steel}}\right).
\]
若锚链接触海床，其悬空段底端竖向张力为零；否则全长悬空。用悬链线垂向落差与各刚体的垂向投影之和构造 \(H(\delta)\)，以 \(H(\delta)=18\) m 求解吃水。

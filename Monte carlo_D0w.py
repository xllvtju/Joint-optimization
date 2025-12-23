# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 15:15:45 2024

@author: 吕晓磊
"""

import numpy as np
import time

# 更新参数定义
R_max = 1  # 最大生产率
lambda_control = 1 / 20  # 控制系统失效率 (泊松分布)
shape_gamma = 0.4  # 伽马过程的形状参数（退化）
beta_0 = 3  # 初始退化率参数,尺度参数
eta = 1.5  # 退化率动态公式参数
gamma_exp = 0.3  # 退化率动态公式参数

# 以下参数将进行敏感性分析
D_0 = 0.02  # 缺陷率初始值
w = 0.68  # 缺陷率参数

u = 2 * np.pi * 10 ** -6  # 缺陷率参数
v = 2  # 缺陷率指数
product_profit = 20  # 每单位产品的利润
a = 0.2  # 生产率函数参数
b = 3    # 生产率函数参数
m, n = 0.2, 0.6  # 退化计算参数
C_PM, C_CM, C_CD, C_PD, C_MR = 200, 500, 300, 100, 100  # 维护成本参数
Lc = 25  # 纠正性维护阈值
q = 20  # 每单位生产率的质量损失成本

# 计算缺陷率函数，添加D_0和w参数接收
def calculate_defect_rate(X_t, D_0, w):
    return D_0 + w * (1 - np.exp(-u * (X_t ** v)))

# 计算退化率 beta 并用于更新 scale_gamma
def calculate_beta(R):
    return beta_0 * eta * (R ** gamma_exp)

# 计算生产率 R
def calculate_production_rate(X_0, IFP):
    return R_max * np.exp(-IFP * a * (X_0 ** b))

# 定义总收益和质量损失成本的计算函数，添加D_0和w参数接收
def calculate_total_revenue(scenario, R, pi, defect_rate_func, T_hat, D_0, w):
    dt = 0.5  # 增大步长以减少迭代次数
    total_revenue, t = 0.0, 0.0
    max_iterations = 1000  # 最大迭代次数防止无限循环
    iterations = 0
    X_t = np.random.gamma(shape_gamma, 1 / beta_0)  # 初始化 X(t)

    while t < T_hat and iterations < max_iterations:
        D_X_t = defect_rate_func(X_t, D_0, w)
        total_revenue += R * pi * (1 - D_X_t) * dt

        # 实时更新 X(t)
        X_t += np.random.gamma(shape=shape_gamma, scale=1 / calculate_beta(R)) * dt
        t += dt
        iterations += 1

    return total_revenue

def calculate_quality_loss_cost(scenario, R, q, defect_rate_func, T_hat, D_0, w):
    dt = 0.3  # 增大步长以减少迭代次数
    total_quality_loss_cost, t = 0.0, 0.0
    max_iterations = 1000  # 最大迭代次数防止无限循环
    iterations = 0
    X_t = np.random.gamma(shape_gamma, 1 / beta_0)  # 初始化 X(t)

    while t < T_hat and iterations < max_iterations:
        D_X_t = defect_rate_func(X_t, D_0, w)
        total_quality_loss_cost += R * D_X_t * q * dt

        # 实时更新 X(t)
        X_t += np.random.gamma(shape=shape_gamma, scale=1 / calculate_beta(R)) * dt
        t += dt
        iterations += 1

    return total_quality_loss_cost

# 修改calculate_maintenance_cost函数定义，使其接收T_PM、T_CM、T_CD参数
def calculate_maintenance_cost(scenario, T_PM, T_CM, T_CD):
    if scenario == 1:
        return C_MR + C_PD
    elif scenario == 2:
        return C_PM + T_PM + C_PD
    elif scenario == 3:
        return C_MR + T_CD + C_CD
    elif scenario == 4:
        return C_PM + max(T_CD, T_PM) + C_CD
    elif scenario == 5:
        return C_CM + T_CM + C_PD

# 模拟控制系统失效时间的函数（指数分布）
def generate_failure_time(lambda_s):
    return np.random.exponential(1 / lambda_s)

# 使用伽马过程模拟退化时间
def simulate_degradation_time(beta, current_degradation, target_degradation, tau):
    total_time = 0
    while current_degradation < target_degradation and total_time < tau:
        increment = np.random.gamma(shape=shape_gamma, scale=1 / beta)
        current_degradation += increment
        total_time += increment
    return total_time if total_time < tau else float('inf')

# 主蒙特卡洛仿真函数，结合系统失效和设备退化过程，添加D_0、w、T_PM、T_CM、T_CD参数接收
def monte_carlo_simulation_cycle(tau, Lp, IFP, q, D_0, w, T_PM, T_CM, T_CD):
    np.random.seed(int(time.time() * 1000) % 2 ** 32)  # 使用时间戳作为随机种子以增加随机性

    lam = lambda_control * tau
    X_0 = np.random.gamma(shape_gamma, 1 / beta_0)

    # 计算生产率 R 并赋值
    R = calculate_production_rate(X_0, IFP)

    beta = calculate_beta(R)
    scale_gamma = 1 / beta

    # 生成控制系统失效的时间 T1（指数分布）
    T1 = generate_failure_time(lambda_control)

    # 使用伽马过程模拟 T2 和 T3
    T2 = simulate_degradation_time(beta, X_0, Lp, tau) if T1 < tau else float('inf')
    T3 = simulate_degradation_time(beta, X_0, Lc, tau)

    # 根据不同情景设置 T_hat 和 E(T)，并限制 T_hat 最大为 100
    if T1 < tau and T2 >= tau and T3 >= tau:
        scenario = 1
        T_hat = min(tau, 100)
        E_T = tau
        X_0 = (m / (m + n)) * X_0  # 更新 X_0
    elif T1 < tau and T2 < tau and T3 >= tau:
        scenario = 2
        T_hat = min(tau, 100)
        E_T = tau + T_PM
        X_0 = 0  # 情景 II 中 X_0 重置为 0
    elif T1 >= tau and T3 < tau:
        scenario = 3
        T_hat = min(T1, 100)
        E_T = T1 + T_CD
        X_0 = (m / (m + n)) * X_0  # 更新 X_0
    elif T1 >= tau and T2 < tau and T3 >= tau:
        scenario = 4
        T_hat = min(T2, 100)
        E_T = T2 + max(T_CD, T_PM)
        X_0 = 0  # 情景 IV 中 X_0 重置为 0
    else:
        scenario = 5
        T_hat = min(T3, 100)
        E_T = T3 + T_CM
        X_0 = 0  # 情景 V 中 X_0 重置为 0

    # 使用计算的生产率 R 进行收益和成本计算，传入D_0和w参数
    total_revenue = calculate_total_revenue(scenario, R, product_profit, calculate_defect_rate, T_hat, D_0, w)
    quality_loss_cost = calculate_quality_loss_cost(scenario, R, q, calculate_defect_rate, T_hat, D_0, w)
    maintenance_cost = calculate_maintenance_cost(scenario, T_PM, T_CM, T_CD)

    net_revenue = total_revenue - (maintenance_cost + quality_loss_cost)

    # 计算单位时间净收益
    net_revenue_per_time = net_revenue / E_T
    return net_revenue_per_time

# 单线程运行蒙特卡洛仿真，添加D_0、w、T_PM、T_CM、T_CD参数传递
def monte_carlo_simulation_with_time_model(t, Lp, IFP, q, D_0, w, T_PM, T_CM, T_CD, num_cycles=300):
    results = [monte_carlo_simulation_cycle(t, Lp, IFP, q, D_0, w, T_PM, T_CM, T_CD) for _ in range(num_cycles)]
    total_net_revenue = sum(results)
    best_net_revenue = max(results)
    return total_net_revenue / num_cycles, best_net_revenue

# 网格搜索算法实现，添加D_0、w、T_PM、T_CM、T_CD参数接收
def grid_search_algorithm(tau_range, IFP_range, Lp_range, D_0, w, T_PM, T_CM, T_CD, num_cycles=150):
    best_solution = None
    best_fitness = float('-inf')

    for tau in tau_range:
        for IFP in IFP_range:
            for Lp in Lp_range:
                _, fitness_value = monte_carlo_simulation_with_time_model(tau, Lp, IFP, q, D_0, w, T_PM, T_CM, T_CD, num_cycles=num_cycles)
                if fitness_value > best_fitness:
                    best_fitness = fitness_value
                    best_solution = (tau, IFP, Lp)

    return best_solution, best_fitness

# 初始的维护时间参数值
T_PM, T_CM, T_CD = 1, 1.5, 1.5  # 维护时间参数

# 设置搜索范围
tau_range = np.arange(20, 40, 1)
IFP_range = np.arange(0, 1.5, 0.1)
Lp_range = np.arange(10, 25, 1)

# 对D_0进行敏感性分析
D_0_values = [0, 0.02, 0.1, 0.2, 0.3]
print("Sensitivity analysis for D_0:")
for D_0_value in D_0_values:
    start_time = time.time()
    best_solution, best_fitness = grid_search_algorithm(tau_range, IFP_range, Lp_range, D_0_value, w, T_PM, T_CM, T_CD)
    end_time = time.time()
    print(f"For D_0 = {D_0_value}, w = {w}:")
    print("Best Solution (tau, IFP, Lp):", best_solution)
    print("Best Net Revenue:", best_fitness)
    print("Execution Time (seconds):", end_time - start_time)
    print("-" * 50)

# 对w进行敏感性分析
w_values = [0.1, 0.2, 0.4, 0.6, 0.8]
print("Sensitivity analysis for w:")
for w_value in w_values:
    start_time = time.time()
    best_solution, best_fitness = grid_search_algorithm(tau_range, IFP_range, Lp_range, D_0, w_value, T_PM, T_CM, T_CD)
    end_time = time.time()
    print(f"For D_0 = {D_0}, w = {w_value}:")
    print("Best Solution (tau, IFP, Lp):", best_solution)
    print("Best Net Revenue:", best_fitness)
    print("Execution Time (seconds):", end_time - start_time)
    print("-" * 50)

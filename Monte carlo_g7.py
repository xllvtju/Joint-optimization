# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 10:44:26 2024

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

D_0 = 0.02  # 缺陷率初始值
w = 0.68  # 缺陷率参数
u = 2 * np.pi * 10**-6  # 缺陷率参数
v = 2  # 缀陷率指数
product_profit = 20  # 每单位产品的利润
a = 0.2  # 生产率函数参数
b = 3    # 生产率函数参数
m, n = 0.2, 0.6  # 退化计算参数
C_PM, C_CM, C_CD, C_PD, C_MR = 200, 500, 300, 100, 100  # 维护成本参数
T_PM, T_CM, T_CD = 1, 1.5, 1.5  # 维护时间参数
Lc = 25  # 纠正性维护阈值
q = 20  # 每单位生产率的质量损失成本

# 计算缺陷率函数
def calculate_defect_rate(X_t):
    return D_0 + w * (1 - np.exp(-u * (X_t ** v)))

# 计算退化率 beta 并用于更新 scale_gamma
def calculate_beta(R):
    return beta_0 * eta * (R ** gamma_exp)

# 计算生产率 R
def calculate_production_rate(X_0, IFP):
    return R_max * np.exp(-IFP * a * (X_0 ** b))

# 定义总收益和质量损失成本的计算函数
def calculate_total_revenue(scenario, R, pi, defect_rate_func, T_hat):
    dt = 1.0  # 增大步长以减少迭代次数
    total_revenue, t = 0.0, 0.0
    max_iterations = 1000  # 最大迭代次数防止无限循环
    iterations = 0
    X_t = np.random.gamma(shape_gamma * t, 1 / beta_0)  # 初始化 X(t)

    while t < T_hat and iterations < max_iterations:
        D_X_t = defect_rate_func(X_t)
        total_revenue += R * pi * (1 - D_X_t) * dt
        X_t += np.random.gamma(shape=shape_gamma * t, scale=1 / calculate_beta(R)) * dt
        t += dt
        iterations += 1

    return total_revenue

def calculate_quality_loss_cost(scenario, R, q, defect_rate_func, T_hat):
    dt = 1.0  # 增大步长以减少迭代次数
    total_quality_loss_cost, t = 0.0, 0.0
    max_iterations = 1000  # 最大迭代次数防止无限循环
    iterations = 0
    X_t = np.random.gamma(shape_gamma * t, 1 / beta_0)  # 初始化 X(t)

    while t < T_hat and iterations < max_iterations:
        D_X_t = defect_rate_func(X_t)
        total_quality_loss_cost += R * D_X_t * q * dt
        X_t += np.random.gamma(shape=shape_gamma * t, scale=1 / calculate_beta(R)) * dt
        t += dt
        iterations += 1

    return total_quality_loss_cost

def calculate_maintenance_cost(scenario):
    if scenario == 1:
        return C_MR + C_PD
    elif scenario == 2:
        return C_PM + C_PD
    elif scenario == 3:
        return C_MR + C_CD
    elif scenario == 4:
        return C_PM + C_CD
    elif scenario == 5:
        return C_CM + C_PD

# 模拟控制系统失效时间的函数（指数分布）
def generate_failure_time(lambda_s):
    return np.random.exponential(1 / lambda_s)

# 使用伽马过程模拟退化时间
def simulate_degradation_time(beta, current_degradation, target_degradation, tau):
    total_time = 0
    max_iterations = 1000  # 限制最大迭代次数
    iterations = 0

    while current_degradation < target_degradation and total_time < tau:
        increment = np.random.gamma(shape=shape_gamma, scale=1 / beta)
        current_degradation += increment
        total_time += increment
        iterations += 1
        if iterations >= max_iterations:
            break

    return total_time if total_time < tau else float('inf')

# 主蒙特卡洛仿真函数，结合系统失效和设备退化过程
def monte_carlo_simulation_cycle(tau, Lp, IFP, q):
    np.random.seed(int(time.time() * 1000) % 2**32)  # 使用时间戳作为随机种子以增加随机性

    X_0 = np.random.gamma(shape_gamma * tau, 1 / beta_0)

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
        T_hat = min(tau, 50)
        E_T = tau
        X_0 = (m / (m + n)) * X_0  # 更新 X_0
    elif T1 < tau and T2 < tau and T3 >= tau:
        scenario = 2
        T_hat = min(tau, 50)
        E_T = tau + T_PM
        X_0 = 0  # 情景 II 中 X_0 重置为 0
    elif T1 >= tau and T3 < tau:
        scenario = 3
        T_hat = min(T1, 50)
        E_T = T1 + T_CD
        X_0 = (m / (m + n)) * X_0  # 更新 X_0
    elif T1 >= tau and T2 < tau and T3 >= tau:
        scenario = 4
        T_hat = min(T2, 50)
        E_T = T2 + max(T_CD, T_PM)
        X_0 = 0  # 情景 IV 中 X_0 重置为 0
    else:
        scenario = 5
        T_hat = min(T3, 50)
        E_T = T3 + T_CM
        X_0 = 0  # 情景 V 中 X_0 重置为 0

    # 使用计算的生产率 R 进行收益和成本计算
    total_revenue = calculate_total_revenue(scenario, R, product_profit, calculate_defect_rate, T_hat)
    quality_loss_cost = calculate_quality_loss_cost(scenario, R, q, calculate_defect_rate, T_hat)
    maintenance_cost = calculate_maintenance_cost(scenario)

    net_revenue = total_revenue - (maintenance_cost + quality_loss_cost)

    # 计算单位时间净收益
    net_revenue_per_time = net_revenue / E_T
    return net_revenue_per_time

# 修复后的 monte_carlo_simulation_with_time_model
def monte_carlo_simulation_with_time_model(t, Lp, IFP, q, num_cycles=2000):
    results = [monte_carlo_simulation_cycle(t, Lp, IFP, q) for _ in range(num_cycles)]
    total_net_revenue = sum(results)
    best_net_revenue = max(results)
    return total_net_revenue / num_cycles, best_net_revenue

# 遗传算法实现
def genetic_algorithm(tau_range, IFP_range, Lp_range, population_size=20, generations=50, mutation_rate=0.1, num_cycles=1000):
    # 初始化种群
    population = np.random.rand(population_size, 3)
    population[:, 0] = population[:, 0] * (tau_range[1] - tau_range[0]) + tau_range[0]  # tau
    population[:, 1] = population[:, 1] * (IFP_range[1] - IFP_range[0]) + IFP_range[0]  # IFP
    population[:, 2] = population[:, 2] * (Lp_range[1] - Lp_range[0]) + Lp_range[0]  # Lp

    best_solution = None
    best_fitness = float('-inf')

    for generation in range(generations):
        fitness_values = np.array([monte_carlo_simulation_with_time_model(tau, Lp, IFP, q, num_cycles)[0] for tau, IFP, Lp in population])

        # 选择最优个体（轮盘赌选择或选择前k个最优）
        best_idx = np.argmax(fitness_values)
        if fitness_values[best_idx] > best_fitness:
            best_fitness = fitness_values[best_idx]
            best_solution = population[best_idx]

        # 选择操作：根据适应度选择父代
        parents = population[np.argsort(fitness_values)[-population_size // 2:]]  # 选择前一半
        new_population = []

        # 交叉操作
        while len(new_population) < population_size:
            p1, p2 = parents[np.random.choice(parents.shape[0], size=2, replace=False)]
            crossover_point = np.random.randint(1, 3)  # 随机选择交叉点
            child1 = np.copy(p1)
            child2 = np.copy(p2)
            child1[crossover_point:], child2[crossover_point:] = p2[crossover_point:], p1[crossover_point:]

            # 变异操作
            if np.random.rand() < mutation_rate:
                child1[np.random.randint(3)] += np.random.normal(0, 0.1)
            if np.random.rand() < mutation_rate:
                child2[np.random.randint(3)] += np.random.normal(0, 0.1)

            # 添加边界约束
            child1[0] = np.clip(child1[0], tau_range[0], tau_range[1])  # tau
            child1[1] = np.clip(child1[1], IFP_range[0], IFP_range[1])  # IFP
            child1[2] = np.clip(child1[2], Lp_range[0], Lp_range[1])  # Lp

            child2[0] = np.clip(child2[0], tau_range[0], tau_range[1])  # tau
            child2[1] = np.clip(child2[1], IFP_range[0], IFP_range[1])  # IFP
            child2[2] = np.clip(child2[2], Lp_range[0], Lp_range[1])  # Lp

            new_population.append(child1)
            new_population.append(child2)

        population = np.array(new_population)[:population_size]

    return best_solution, best_fitness

# 设置搜索范围
tau_range = np.array([10, 40])
IFP_range = np.array([0, 2])
Lp_range = np.array([10, 25])

# 运行遗传算法并测量运行时间
start_time = time.time()
best_solution, best_fitness = genetic_algorithm(tau_range, IFP_range, Lp_range)
end_time = time.time()

# 输出最佳解及运行时间
print("Best Solution (tau, IFP, Lp):", best_solution)
print("Best Net Revenue:", best_fitness)
print("Execution Time (seconds):", end_time - start_time)

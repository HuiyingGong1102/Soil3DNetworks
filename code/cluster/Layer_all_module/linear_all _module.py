import numpy as np
import pandas as pd
import torch
import torch.optim as optim
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from torch.distributions import MultivariateNormal
from xlwt.ExcelMagic import ptgInt

# from scipy.stats import iqr
import json
# import gc
import warnings
from tqdm import tqdm  # 导入tqdm库
# 检查可用的GPU数量
num_gpus = torch.cuda.device_count()
print(f"Number of GPUs available: {num_gpus}")
warnings.filterwarnings("ignore")
# 检查是否有可用的GPU
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)

data_df = pd.read_csv('log_data_module.csv', index_col=0)

# 将Pandas数据框转换为Numpy数组
data_np = data_df.to_numpy().T

# 将Numpy数组转换为PyTorch张量并移动到GPU
data = torch.tensor(data_np, dtype=torch.float32).to(device)
Time_df = pd.read_csv('log_time_module.csv', index_col=0)
Time = torch.tensor(Time_df.to_numpy(), dtype=torch.float32).to(device)


# def dmvnorm(x, mean, sigma, log=True):
#     # 使用 MultivariateNormal 来计算概率密度
#     dist = MultivariateNormal(mean, sigma)
#
#     # 返回对数概率密度或概率密度
#     if log:
#         return dist.log_prob(x)  # 对数概率密度
#     else:
#         return dist.pdf(x)  # 概率密度


def mahalanobis(x, center, cov):
    x_cen = x - center
    sol = torch.linalg.solve(cov, x_cen.T).T
    dist = torch.sum(x_cen * sol, dim=1)
    return dist


def dmvnorm(x, mean, sigma, log=True):
    distval = mahalanobis(x, mean, sigma)
    # 用 Cholesky 分解计算 logdet，更稳定
    L = torch.linalg.cholesky(sigma)
    logdet = 2 * torch.sum(torch.log(torch.diag(L)))

    d = x.shape[1]
    log2pi = torch.log(torch.tensor(2 * np.pi, device=x.device, dtype=torch.float32))

    logretval = -0.5 * (d * log2pi + logdet + distval)

    if log:
        return logretval
    else:
        return torch.exp(logretval)

def linear_equation(x, linear_par):
    result = linear_par[:, 0][:, None] + linear_par[:, 1][:, None] * x
    return result


def linear_equation_base(x, y):
    x = np.array(x, dtype=float).reshape(-1, 1)
    y = np.array(y, dtype=float)

    model = LinearRegression().fit(x, y)
    intercept = model.intercept_
    slope = model.coef_[0]

    return intercept, slope


def logsumexp(v):
    vm = torch.max(v)
    return torch.log(torch.sum(torch.exp(v - vm))) + vm


# def get_SAD1_covmatrix(par, n):
#     device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#
#     # 转到目标 device
#     par = torch.as_tensor(par, dtype=torch.float32, device=device)
#
#     phi = par[0]
#     gamma = par[1]
#
#     # 创建索引差矩阵
#     idx = torch.arange(n, dtype=torch.float32, device=device)
#     diff_matrix = torch.abs(idx[:, None] - idx[None, :])
#
#     # 计算 AR(1) 协方差
#     sigma = (phi ** diff_matrix) * (gamma ** 2)
#
#     # 强制对称化
#     sigma = (sigma + sigma.T) / 2
#
#     # 处理 NaN / Inf
#     sigma = torch.nan_to_num(sigma, nan=0.0, posinf=0.0, neginf=0.0)
#
#     # 加稳定项
#     eps = 1e-6
#     sigma = sigma + eps * torch.eye(n, dtype=torch.float32, device=device)
#
#     return sigma

# 辅助函数：计算协方差矩阵
def get_SAD1_covmatrix(par, n):
    phi = par[0]
    gamma = par[1]

    # 创建索引矩阵
    indices = torch.arange(1, n + 1, dtype=torch.float32).to(device)

    # 计算对角线元素
    diag_elements = (1 - torch.pow(phi, 2 * indices)) / (1 - torch.pow(phi, 2))

    # 计算非对角线元素
    sigma = torch.zeros((n, n), dtype=torch.float32).to(device)
    for i in range(n):
        sigma[i, i:] = torch.pow(phi, torch.arange(n - i, dtype=torch.float32).to(device)) * diag_elements[i]

    # 使矩阵对称
    sigma = sigma + sigma.T - torch.diag(torch.diag(sigma))

    # 缩放矩阵
    sigma *= torch.pow(gamma, 2)
    eps = 1e-6
    sigma = sigma + eps * torch.eye(n, dtype=torch.float32, device=device)
    return sigma

# def get_SAD1_covmatrix(par, n):
#     eps = 1e-6
#     phi = par[0]
#     gamma = par[1]
#     indices = torch.arange(1, n + 1, dtype=torch.float32).to(device)
#
#     # 计算对角元素
#     diag_elements = (1 - torch.pow(phi, 2 * indices)) / (1 - torch.pow(phi, 2))
#
#     # 只创建对角矩阵，不包含非对角元素
#     sigma = torch.diag(diag_elements) * torch.pow(gamma, 2)
#
#     return sigma + eps * torch.eye(n, dtype=torch.float32, device=device)




def get_par_int(X, k, times1, n1):
    n, d = X.shape
    X1 = X.cpu()

    cov_int = [0.1, 0.1]

    # x_t = X[:, 1:]
    # x_tm1 = X[:, :-1]
    #
    # num = torch.sum(x_t * x_tm1)
    # den = torch.sum(x_tm1 * x_tm1) + 1e-12
    # phi_hat = num / den
    # phi_hat = torch.clamp(phi_hat, -0.9999, 0.9999)
    #
    # resid = x_t - phi_hat * x_tm1
    # gamma2_hat = torch.sum(resid * resid) / (n * (d - 1))
    # gamma_hat = torch.sqrt(gamma2_hat + 1e-12)
    #
    # cov_int = [float(phi_hat), float(gamma_hat)]

    print(cov_int)
    init_cluster = KMeans(n_clusters=k,
                          init='k-means++',
                          n_init=10
                          ).fit(X.cpu())
    labels = init_cluster.labels_

    prob = np.bincount(labels) / n
    times1 = times1.cpu().numpy()

    fit1 = np.array([linear_equation_base(times1, init_cluster.cluster_centers_[c, :n1]) for c in range(k)]).T
    return_obj = {
        'initial_cov_params': cov_int,
        'initial_mu_params': fit1.flatten(),
        'initial_probibality': prob
    }

    return return_obj


# 主函数
def Q_function(par, prob_log, omega_log, X, k, n1, times1):
    n = X.shape[0]
    n1, k = map(int, [n1, k])

    X1 = X

    par_mu = par[2:]

    cov1 = get_SAD1_covmatrix(par[:2], n1)

    mu1_mat = torch.column_stack((par_mu[:k], par_mu[k:2 * k]))

    mu1 = linear_equation(times1, mu1_mat)

    mvn_log1 = torch.stack([dmvnorm(X1, mu1[i], cov1, True) for i in range(k)], dim=1)
    mvn_log = mvn_log1
    tmp = prob_log + mvn_log - omega_log
    Q = -torch.sum(tmp * torch.exp(omega_log))

    return Q


def fun_clu(data, k, Time, n1, initial_pars=None,iter_max=10000):

    d, n = data.shape
    delta_per = 10000
    rel_impr = 10000
    iter = 0

    X1 = data
    times1 = Time

    if initial_pars is None:
        par_int = get_par_int(data, k, times1, n1)
        prob_int = np.ones(k) / k
        # prob_int = par_int['initial_probibality']
        cov_int = par_int['initial_cov_params']
        mu_int = par_int['initial_mu_params']
        initial_pars = np.hstack((cov_int, mu_int))
    else:
        prob_int = torch.tensor(initial_pars[:k]).to(device)
        initial_pars = initial_pars[k:]

    par = torch.tensor(initial_pars, dtype=torch.float32).to(device)
    par.requires_grad = True
    prob_logi = torch.log(torch.tensor(prob_int, dtype=torch.float32).to(device))
    # prob_logi = torch.log(torch.tensor(prob_int, dtype=torch.float32))

    initial_lr = 1e-8
    optimizer = optim.AdamW([par], lr=initial_lr)

    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=3)

    while abs(delta_per) > 0.001 and abs(rel_impr) > 0.0001 and iter < iter_max:
        par_hat = par
        par_mui = par[2:]
        cov1i = get_SAD1_covmatrix(par[:2], n1).to(device).detach()

        mu1_mati = torch.column_stack((par_mui[:k], par_mui[k:2 * k])).to(device).detach()

        mu1i = linear_equation(times1, mu1_mati).to(device).detach()  # 需要的参数是按行排列的

        mvn_log1i = torch.stack([dmvnorm(X1, mu1i[i], cov1i, True) for i in range(k)], dim=1).to(device).detach()

        mvn_logi = mvn_log1i
        mvni = mvn_logi + prob_logi
        omega_logi = mvni - torch.logsumexp(mvni, dim=1, keepdim=True)
        omegai = torch.exp(omega_logi)

        # omega_logi = torch.stack([mvni[i] - torch.logsumexp(mvni[i], dim=0) for i in range(d)])

        LL_mem = Q_function(par, prob_logi, omega_logi, data, k, n1, times1)
        print(LL_mem)

        nk = omegai.sum(dim=0)
        wk = nk.clamp_min(1e-12)
        pi = wk.pow(1.0 / max(2, 1e-8))
        pi = pi / pi.sum()

        # 平滑
        alpha = 1 / (2 * k)
        prob = pi + alpha
        prob = prob / prob.sum()

        prob_logi = torch.log(prob)

        # prob_expi = torch.stack([logsumexp(omega_logi[:, i]) for i in range(omega_logi.size(1))])
        # prob_logi = prob_expi - torch.log(torch.tensor(n1, dtype=torch.float32))
        # optimizer = optim.AdamW([par],lr=1e-10)

        for _ in tqdm(range(25), desc="Training Progress"):  # 为训练循环添加进度条
            optimizer.zero_grad()
            Q = Q_function(par, prob_logi, omega_logi, data, k, n1, times1)
            Q.backward()
            torch.nn.utils.clip_grad_norm_(par, max_norm=1.0)
            optimizer.step()

        num = torch.norm(par_hat - par)  # Frobenius/L2
        den = torch.norm(par_hat) + 1e-12

        par_hat = par
        par = par_hat

        LL_next = Q_function(par, prob_logi, omega_logi, data, k, n1, times1)
        # epsilon = (LL_next - LL_mem)/LL_next
        delta = LL_next - LL_mem
        delta_per = delta / (d + 1e-12)  # 每样本改进
        rel_impr = delta / (abs(LL_mem) + 1e-12)  # 相对改进
        LL_mem = LL_next
        scheduler.step(LL_next)
        iter += 1

        print(f"Iter: {iter}, Log-Likelihood: {LL_next}")
        final_assignments = torch.argmax(omegai, dim=1).cpu().numpy()
        print(f"各聚类的样本数量: {np.bincount(final_assignments, minlength=k)}")

    prob_log = prob_logi.detach().cpu().numpy()
    omega_logi_cpu = omega_logi.cpu()
    max_values, max_indices = torch.max(omega_logi_cpu, dim=1)
    max_omega_logi = max_indices.detach().numpy()
    par = par.detach().cpu().numpy()
    LL_next = LL_next.detach().cpu().numpy()
    BIC = 2 * (LL_next) + np.log(d) * (len(par) + k - 1)
    print(par[0:2])
    return {
        'cluster_number': k,
        'log-likelihood': LL_next,
        'BIC': BIC,
        'par': par,
        'prob_log': prob_log,
        'max_omega_logi': max_omega_logi
    }


import json
import os


# def run_analysis(k, data):
#     """
#     执行单个 k 值的聚类分析，并保存结果到 JSON 文件。
#     """
#     print(f"开始运行 k={k} 的分析...")
#     result = fun_clu(data, k, Time, n1=90)
#     result['par'] = result['par'].tolist()
#     result['log-likelihood'] = result['log-likelihood'].tolist()
#     result['prob_log'] = result['prob_log'].tolist()
#     result['max_omega_logi'] = result['max_omega_logi'].tolist()
#
#     # 动态生成保存路径
#     save_path = f'./json/k{k}.json'
#     # save_path = f'./json/k3-{i}.json'
#
#     # 确保路径存在
#     os.makedirs(os.path.dirname(save_path), exist_ok=True)
#
#     # 保存结果到 JSON 文件
#     with open(save_path, 'w') as f:
#         json.dump(result, f)
#
#     print(f"k={k} 的分析完成，结果保存在 {save_path}")


# k_values = list(range(155,185,5))
# k_values = [60]
#  # 设置需要分析的 k 值范围
# for k in k_values:
#     run_analysis(k, data)


for k in list(range(35,38,1)):
    for i in list(range(31, 51, 1)):
        print(f"开始运行 k={k} 的分析...")
        result = fun_clu(data, k, Time, n1=90)
        result['par'] = result['par'].tolist()
        result['log-likelihood'] = result['log-likelihood'].tolist()
        result['prob_log'] = result['prob_log'].tolist()
        result['max_omega_logi'] = result['max_omega_logi'].tolist()

        # 动态生成保存路径
        # save_path = f'./json/k{k}.json'
        save_path = f'./json_module-0.1/k{k}-{i}.json'

        # 确保路径存在
        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        # 保存结果到 JSON 文件
        with open(save_path, 'w') as f:
            json.dump(result, f)

        print(f"k={k} 的分析完成，结果保存在 {save_path}")
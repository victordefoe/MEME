
# This is a project to predict the tendancy of WuHan 19 nCoVirus
# The 19_nCoV source data are from the project:
# https://github.com/globalcitizen/2019-wuhan-coronavirus-data

import pandas as pd
import os
from os.path import join as opj
import re
import time
import datetime
import json
import numpy as np
import matplotlib.pyplot as plt
import geatpy as ea
from utils import residual_square
import utils
import math

def rough_data_clean():
    time_pt = re.compile(r'[0-9]') 
    dp = './19_nCoV_data/data-sources/dxy/data/'
    namelist = os.listdir(dp)
    namelist = [name for name in namelist if name[-3:] == 'csv']
    for n in namelist:
        data_org = opj(dp, n)
        save_path = opj('./19_nCoV_clean_data/rough',n[:15]+'.csv')
        a = pd.read_csv(data_org, header=None)
        # time = time_pt.findall(list(a.iloc[1])[0])
        # time_day = '-'.join([''.join(time[:4]),''.join(time[4:6]), ''.join(time[6:8])])
        # time_sec = ':'.join([''.join(time[8:10]), ''.join(time[10:12]), ''.join(time[12:])])
        # time = time_day + ' ' + time_sec
        time = n[:15]
        # time = ' '.join(['-'.join([time[0][:4], time[0][4:6], time[0][6:]]), 
        #         ':'.join([time[1][:2], time[1][2:4], time[1][4:]])])
        time = datetime.datetime.strptime(time, "%Y%m%d-%H%M%S")


        a = a.iloc[3:-1].reset_index(drop=True)

        pp = [[], [], []]
        for each in a[0]:
            pp[0].append(each.split('|')[0])
            pp[1].append(int(each.split('|')[1]))
            pp[2].append(int(each.split('|')[2]))
        pp[0].append('China Total')
        pp[1].append(sum(pp[1]))
        pp[2].append(sum(pp[2]))
        places = pd.Series(pp[0])
        confs = pd.Series(pp[1])
        deaths = pd.Series(pp[2])
        new_df = pd.DataFrame(columns=['place', 'confirmed','deaths'])
        new_df['place'] = places
        new_df['confirmed'] = confs
        new_df['deaths'] = deaths

        new_df.to_csv(save_path, header=True, index=False)


        ## save all the death data
        
def detail_data_clean(prov_shortname,cityname):
    # para: prov_shortname: shortname of the province you are predicting
    # para: cityname: shortname of the city you are predicting
    dp = './19_nCoV_data/data-sources/dxy/data/'
    namelist = os.listdir(dp)
    namelist = [name for name in namelist if name[-4:] == 'json']
    times, confirmedCount, suspectedCount, curedCount, deadCount = [], [], [], [], []

    for n in namelist:
        time = n[:15]
        time = datetime.datetime.strptime(time, "%Y%m%d-%H%M%S")
        times.append(time)
        with open(opj(dp, n), 'r', encoding='utf-8') as f:
            data = json.load(f)
        for eachprov in data:
            if eachprov['provinceShortName'] != prov_shortname:
                continue
            for eachcity in eachprov['cities']:
                if eachcity['cityName'] != cityname:
                    continue
                # Place.append(eachprov['cityName'])
                confirmedCount.append(eachcity['confirmedCount'])
                suspectedCount.append(eachcity['suspectedCount'])
                curedCount.append(eachcity['curedCount'])
                deadCount.append(eachcity['deadCount'])
    city_data = pd.DataFrame(columns=[
        'time', 'confirmedCount', 'suspectedCount', 'curedCount','deadCount'])
    satistics = [times, confirmedCount, suspectedCount, curedCount, deadCount]
    
    for idx, eachcol in enumerate(list(city_data.columns)):
        city_data[eachcol] = pd.Series(satistics[idx])
    return city_data


    
class E_SEIR(object):
    """
    build up an E-SEIR model for 19_nCoV plague, help with projects following:
    https://github.com/XuelongSun/Dynamic-Model-of-Infectious-Diseases
    """
    def __init__(self, T=12):
        # population
        self.N = 14186500  # data from the website
        # simuation Time / Day
        self.T = T
        # susceptiable ratio
        self.s = np.zeros([self.T])
        # exposed ratio
        self.e = np.zeros([self.T])
        # infective ratio
        self.i = np.zeros([self.T])
        # remove ratio
        self.r = np.zeros([self.T])

        # contact rate
        self.lamda = 0  # = 感染者日均接触人数 * 接触传染概率 / 总人数
        self.omiga = 0.1 # effect the lower rate of the lamda
        self.d = 20 # deviation of time
        # isolate rate 患病的人中被隔离的比率 (zf措施)
        self.beta = 0.8
        # remove rate # 死亡率+治愈率
        self.gamma = 0.3
        # exposed period # 平均潜伏期的倒数
        self.sigma = 1 / 7

        # initial infective people
        self.i[0] = 1.0 / self.N
        self.s[0] = (1e7 + 4e6) / self.N
        self.e[0] = 10.0 / self.N
    
    def deduce(self):
        for t in range(self.T-1):
            self.lamda = 1 - 1 / (1 + math.exp(-1 * self.omiga * (t-self.d)))
            self.s[t + 1] = self.s[t] \
                - self.lamda * self.s[t] * (self.i[t] * (1- self.beta) + self.e[t])
            self.e[t + 1] = self.e[t] \
                + self.lamda * self.s[t] * (self.i[t] * (1-self.beta) + self.e[t]) \
                            - self.sigma * self.e[t]
            self.i[t + 1] = self.i[t] + self.sigma * self.e[t] - self.gamma * self.i[t]
            self.r[t + 1] = self.r[t] + self.gamma * self.i[t]

    def draw_curves(self, start_day=0, end_day=1000):
        self.s *= self.N
        self.e *= self.N
        self.i *= self.N
        self.r *= self.N

        if end_day > self.T:
            end_day = self.T
    
        fig, ax = plt.subplots(figsize=(10,6))
        # ax.plot(self.s[start_day:end_day], c='b', lw=2, label='S')
        ax.plot(self.e[start_day:end_day], c='orange', lw=2, label='E')
        ax.plot(self.i[start_day:end_day], c='r', lw=2, label='I')
        ax.plot(self.r[start_day:end_day], c='g', lw=2, label='R')
        ax.set_xlabel('Day',fontsize=15)
        ax.set_ylabel('Number of people', fontsize=15)
        ax.grid(1)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.legend()
        plt.show()





def train_process(model, gt):
    def aim(Phen):
        # ld = Phen[:, [0]] # 取出第1列，得到所有个体的第1个自变量
        # gm = Phen[:, [1]] # 取出第2列，得到所有个体的第2个自变量
        fits = np.zeros((len(Phen),1))  # 所有个体的健康程度评估
        # print(Phen)
        for idx in range(len(Phen)):
            model.omiga = Phen[idx][1] 
            model.gamma = Phen[idx][0]
            # model.beta = Phen[idx][2]
            model.d = Phen[idx][2]
             
            model.deduce()
            # print(np.array(model.i * model.N))
            exam_day = int(Phen[idx][3])

            model_preds_i = model.i[exam_day:exam_day+len(gt)] * model.N
            gt_i = np.array(gt['confirmedCount'])
            model_preds_r = model.r[exam_day:exam_day+len(gt)] * model.N
            gt_r = np.array(gt['deadCount']+gt['curedCount'])

            infective_fitness = residual_square(model_preds_i,  gt_i)
            recov_fitness = residual_square(model_preds_r,  gt_r)
            
            fits[idx] = infective_fitness + recov_fitness
        # print(model_preds_i, gt_i)
        # print(fits)
        return fits
          

    # settings
    omiga = [0, 1]                   # 自变量范围
    gamma = [0, 1]
    # beta = [0, 1]
    d = [0, 100]
    exam_d = [10, 60]

    b1 = [0, 0]                  # 自变量边界, 1 表示包含边界， 0 表示不包含边界
    b2 = [0, 0]
    # b3 = [0, 0]
    b4 = [0, 0]
    b5 = [1, 1]
    varTypes = np.array([0, 0, 0, 0])       # 自变量的类型，0表示连续，1表示离散
    Encoding = 'BG'                # 'BG'表示采用二进制/格雷编码
    codes = [1, 1, 1, 1]                    # 变量的编码方式，2个变量均使用格雷编码
    precisions =[6, 6, 6, 6]                # 变量的编码精度
    scales = [0, 0, 0, 0]                   # 采用算术刻度
    ranges=np.vstack([gamma, omiga, d, exam_d]).T       # 生成自变量的范围矩阵
    borders=np.vstack([b2, b1, b4, b5]).T      # 生成自变量的边界矩阵

    
    # params of GA
    NIND = 400                     # 种群个体数目
    MAXGEN = 100                 # 最大遗传代数
    maxormins = [1]   # 最小化目标函数， 元素为-1表示最大化目标函数
    FieldD = ea.crtfld(Encoding,varTypes,
                        ranges, borders,precisions,codes,scales) # 调用函数创建区域描述器
    Lind = int(np.sum(FieldD[0, :]))          # 计算编码后的染色体长度
    obj_trace = np.zeros((MAXGEN, 2))         # 定义目标函数值记录器
    var_trace = np.zeros((MAXGEN, Lind))      # 定义染色体记录器，记录每一代最优个体的染色体


    start_time = time.time()                             # 开始计时
    Chrom = ea.crtbp(NIND, Lind)                         # 生成种群染色体矩阵
    variable = ea.bs2real(Chrom, FieldD)                 # 对初始种群进行解码
    ObjV = aim(variable)          # 计算初始种群个体的目标函数值
    best_ind = np.argmin(ObjV)                         # 计算当代最优个体的序号


    # 开始进化
    for gen in range(MAXGEN):
        
        FitnV = ea.ranking(maxormins * ObjV)               # 根据目标函数大小分配适应度值(由于遵循目标最小化约定，因此最大化问题要对目标函数值乘上-1)
        SelCh=Chrom[ea.selecting('rws', FitnV, NIND-1), :] # 选择，采用'rws'轮盘赌选择
        SelCh=ea.recombin('xovsp', SelCh, 0.7)           # 重组(采用两点交叉方式，交叉概率为0.7)
        SelCh=ea.mutbin(Encoding, SelCh)                 # 二进制种群变异
        # 把父代精英个体与子代合并
        Chrom = np.vstack([Chrom[best_ind, :], SelCh])
        variables = ea.bs2real(Chrom, FieldD)             # 对育种种群进行解码(二进制转十进制)
        ObjV = aim(variables)                             # 求育种个体的目标函数值
        # 记录
        best_ind = np.argmin(ObjV)                       # 计算当代最优个体的序号
        obj_trace[gen, 0] = np.sum(ObjV) / NIND          # 记录当代种群的目标函数均值
        obj_trace[gen, 1] = ObjV[best_ind]               # 记录当代种群最优个体目标函数值
        var_trace[gen, :] = Chrom[best_ind, :]           # 记录当代种群最优个体的变量值
    # 进化完成
    end_time = time.time() # 结束计时

    # 绘制图像
    ea.trcplot(obj_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])

    best_gen = np.argmin(obj_trace[:, [1]])
    print('最优解的目标函数值: ', obj_trace[best_gen, 1])
    opt_variable = ea.bs2real(var_trace[[best_gen], :], FieldD)
    print('最优解的决策变量为：')
    for i in range(opt_variable.shape[1]):
        print(opt_variable[0, i])
    print('用时: ', end_time - start_time, '秒')
    
    model.omiga = opt_variable[0, 1]
    model.gamma = opt_variable[0, 0]
    # model.beta = opt_variable[0, 2]
    model.d = opt_variable[0, 2]
    model.deduce()
    model.draw_curves()

    return opt_variable[0]

def get_evaluation(model, fit_data, pre_gt, T=150):
    assert T >= 2 * len(fit_data+pre_gt), 'Should set a more long period'

    opt_vars = train_process(model, fit_data) # [gamma, omiga, d, exam_day]


    test_model = E_SEIR(T=T)
    test_model.gamma = opt_vars[0]
    test_model.omiga = opt_vars[1]
    test_model.d = opt_vars[2]
    exam_day = opt_vars[3]
    exam_day = int(exam_day)+len(fit_data)


    test_model.deduce()
    test_model.draw_curves()

    pre_model_i = test_model.i[exam_day: exam_day+len(pre_gt)]
    pre_model_r = test_model.r[exam_day: exam_day+len(pre_gt)]
    pre_gt_i = np.array(pre_gt['confirmedCount'])
    pre_gt_r = np.array(pre_gt['deadCount']+ pre_gt['curedCount'])

    i_error = utils.related_error_rate(pre_model_i, pre_gt_i)
    r_error = utils.related_error_rate(pre_model_r, pre_gt_r)

    print('误差率---  感染者误差率:{}, 移除者误差率:{}'.format(i_error, r_error))
    print(pre_gt_i, pre_model_i)





    



def get_R0(t, Yt, Tg, Ti):
    # t: 初始病发至今，某一时刻
    # Yt：t时刻对应病患数
    # Tg: 潜伏期时长
    # Ti: 感染期时长
    lamda = math.log(Yt) / t
    r0 = (1 + Tg * lamda) * (1 + Ti * lamda)
    return r0

if __name__ == '__main__':
    data = detail_data_clean('湖北','武汉')
    # only save the data on day
    data['time'] = data['time'].dt.date
    data_day = data.groupby(by=['time']).head(1).reset_index(drop=True)

    # confirmed count = accumulated - removed
    data_day['confirmedCount'] = data_day['confirmedCount'] \
        - (data_day['curedCount']+data_day['deadCount'])

    # 提前n天划分数据集
    n = 7
    train_day = data_day.iloc[:-n]
    test_day = data_day.iloc[-n:]
    print(train_day)

    model = E_SEIR(T=150)
    print('划分数据集 训练天数和测试天数分别为：',train_day.shape[0], test_day.shape[0])
    get_evaluation(model, train_day, test_day)
    # model.deduce()
    # fitness = residual_square(np.array(model.i), np.array(data_day['confirmedCount']))
    # model.draw_curves()
    
    R0 = get_R0(40, 572, 7, 1)
    print('基本再生指数：', R0)
    
    
    

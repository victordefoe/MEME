
### repo closed！
>With the policy of China, co-virus will be less care and discussions, the repo closed.    
（疫情已成为过去式，关于最初的感染预测模型也已经应该被封存在历史的长河中。但愿下一次风暴来临的时候，我们能够记住曾经发生过的和战斗过的一切。）

--- 
## MEME: Mathematical Epidemiological Model using Eovolution algorithmn for the COVID-19 virus
Establishing an Epidemiological Model Prediction of COVID-19 Novel Coronavirus Based on Mathematical Model and Optimization Method

基于数学模型和优化方法建立新型冠状病毒的传染病学模型预测

Patient population data obtained from this project   
(患病人数资料由这个项目获取)   
https://github.com/globalcitizen/2019-wuhan-coronavirus-data


You can git clone the above link to the 19_nCoV_data folder  
(可以git clone 上述链接 到19_nCoV_data文件夹下)
```bash
git clone https://github.com/globalcitizen/2019-wuhan-coronavirus-data ./19_nCoV_data/
```

usage (使用方法)   
```bash 
python nCov_19.py
```

### requirments referenece
Some open source packages for the Python platform  
（ 用的Python 平台一些开源包）
```bash     
pandas>=1.2.4
geatpy>=2.7.0  
matplotlib>=3.3.4 
numpy>=1.23.0
json5>=0.9.5
```

# Theory
Improved on the basis of the SEIR model

Based on SEIR, the following considerations are made:
* Considering that the virus is equally contagious during the incubation period
* Since the confirmed patients were better isolated in the hospital, the infection rate of this part of the patients is adjusted considering the isolation coefficient
* Calculated using the average incubation period of 7 days
* Using the s-type function with parameters to fit the infection rate decline law to simulate increasingly perfect epidemic prevention
* Optimization of model parameters through evolutionary calculations. The results show that the current trend (February 15, 2020) is consistent with the trend over the last 36 days in the contagion model.

# Report
The Chinese version of the very detailed report with technical details has been uploaded, see 'docs/report.pdf'


# 理论  
基于SEIR 模型
改进

基于SEIR基础上，做了如下考虑：
* 考虑病毒在潜伏期也具有相同传染性
* 考虑确诊患者在医院进行了较好的隔离，这部分患者比率传染性配比了隔离系数
* 取潜伏期平均日7日计算
* 使用带参数的s型函数拟合传染率的下降规律，来模拟日益完善的防疫力量
* 用演化计算优化拟合参数，结果显示目前趋势(2020年2月15日)相当于模型发病前36天的趋势


# 论文
技术细节和相关论文说明报告已经更新，请看'docs/report.pdf'



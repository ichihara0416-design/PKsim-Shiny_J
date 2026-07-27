# PKsim Shiny
# Pharmacokinetic simulation application based on a linear
# one-compartment model
#
# Copyright (c) 2026 Yoshinori Ichihara
# Licensed under the MIT License.
# See the LICENSE file for details.
#
# This software may be used, modified, and redistributed for
# educational, research, and other non-clinical purposes.
#
# IMPORTANT:
# This software is not a validated medical device and is not intended
# for clinical decision-making, patient care, diagnosis, treatment,
# or individual dose adjustment.
# Do not use its outputs for clinical purposes.
#
# 本ソフトウェアは、教育、研究、その他の非臨床目的で
# 利用、改変および再配布できます。
#
# 重要：
# 本ソフトウェアは承認または検証された医療機器ではありません。
# 実際の患者の診療、診断、治療、投与設計または臨床判断に、
# 本ソフトウェアの出力を使用しないでください。


## Run Appをクリックして、アプリを起動してください。
## 以下アプリ起動用コード：改変の必要はありません。


# 必要なパッケージのインストール
required_packages <- c("shiny", "ggplot2", "markdown", "showtext", "sysfonts", "curl")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(new_packages)) {
  install.packages(new_packages)
}

library(shiny)
library(ggplot2)
library(markdown)
library(showtext)
library(sysfonts)

font_add_google("Noto Sans JP", "jp")
showtext_auto()

# ----------------------------
# 関数
# ----------------------------

# eGFRとCcrの計算関数
calculate_eGFR <- function(scr, age, weight, height, sex) {
  bsa <- 0.007184 * (weight^0.425) * (height^0.725)
  if (sex == "男性") {
    return(175 * (scr ^ -1.154) * (age ^ -0.203) * (bsa / 1.73))
  } else {
    return(175 * (scr ^ -1.154) * (age ^ -0.203) * 0.742 * (bsa / 1.73))
  }
}

calculate_Ccr <- function(scr, age, weight, sex) {
  if (sex == "男性") {
    return(((140 - age) * weight) / (72 * scr))
  } else {
    return((((140 - age) * weight) / (72 * scr)) * 0.85)
  }
}

# t1/2からkeへの変換
t_half_to_ke <- function(t_half) {
  log(2) / t_half
}

# 24時間ごとのAUC
calc_auc24 <- function(Time, Cp) {
  dt <- diff(Time)[1]
  points_per_24h <- round(24 / dt)
  n_block <- floor(length(Time) / points_per_24h)
  
  if (n_block < 1) return(numeric(0))
  
  auc24 <- numeric(n_block)
  for (i in seq_len(n_block)) {
    start_idx <- (i - 1) * points_per_24h + 1
    end_idx <- i * points_per_24h
    auc24[i] <- sum(Cp[start_idx:end_idx]) * dt
  }
  auc24
}

# lower bound crossing time
calc_crossing_times <- function(Time, Cp, lower) {
  idx <- which(diff(sign(Cp - lower)) != 0)
  if (length(idx) == 0) return(numeric(0))
  sort(Time[idx])
}

# 単回 IV bolus
simulate_iv_bolus_single <- function(dose, Vd, ke, Time) {
  Cp <- dose / Vd * exp(-ke * Time)
  AUC <- dose / (ke * Vd)
  list(Cp = Cp, AUC = AUC)
}

# 単回 経口
simulate_oral_single <- function(dose, Vd, ke, ka, F, Time) {
  Cp <- F * dose / Vd * (ka / (ka - ke)) * (exp(-ke * Time) - exp(-ka * Time))
  AUC <- F * dose / (ke * Vd)
  list(Cp = Cp, AUC = AUC)
}

# 点滴静注
simulate_iv_infusion <- function(rate, Vd, ke, infusion_time, Time) {
  Cp <- ifelse(
    Time <= infusion_time,
    rate / (ke * Vd) * (1 - exp(-ke * Time)),
    rate / (ke * Vd) * (1 - exp(-ke * infusion_time)) * exp(-ke * (Time - infusion_time))
  )
  AUC <- sum(Cp) * diff(Time)[1]
  list(Cp = Cp, AUC = AUC)
}

# IV繰り返し
simulate_iv_bolus_repeated <- function(dose, Vd, ke, tau, n, Time) {
  Cp <- rep(0, length(Time))
  AUC_list <- numeric(n)
  
  for (i in seq_len(n)) {
    start_time <- (i - 1) * tau
    Cp_dose <- dose / Vd * exp(-ke * (Time - start_time)) * (Time >= start_time)
    Cp <- Cp + Cp_dose
    AUC_list[i] <- dose / (ke * Vd)
  }
  
  list(
    Cp = Cp,
    AUC_list = AUC_list,
    totalAUC = sum(AUC_list)
  )
}

# IV繰り返し（途中で投与量変更）
simulate_iv_bolus_repeated_changed <- function(initial_dose, changed_dose, initial_dose_num, Vd, ke, tau, n, Time) {
  Cp <- rep(0, length(Time))
  AUC_list <- numeric(n)
  
  for (i in seq_len(n)) {
    this_dose <- ifelse(i <= initial_dose_num, initial_dose, changed_dose)
    start_time <- (i - 1) * tau
    Cp_dose <- this_dose / Vd * exp(-ke * (Time - start_time)) * (Time >= start_time)
    Cp <- Cp + Cp_dose
    AUC_list[i] <- this_dose / (ke * Vd)
  }
  
  list(
    Cp = Cp,
    AUC_list = AUC_list,
    totalAUC = sum(AUC_list)
  )
}

# 経口繰り返し
simulate_oral_repeated <- function(dose, Vd, ke, ka, F, tau, n, Time) {
  Cp <- rep(0, length(Time))
  AUC_list <- numeric(n)
  
  for (i in seq_len(n)) {
    start_time <- (i - 1) * tau
    Cp_dose <- F * dose / Vd * (ka / (ka - ke)) *
      (exp(-ke * (Time - start_time)) - exp(-ka * (Time - start_time))) *
      (Time >= start_time)
    Cp <- Cp + Cp_dose
    AUC_list[i] <- F * dose / (ke * Vd)
  }
  
  list(
    Cp = Cp,
    AUC_list = AUC_list,
    totalAUC = sum(AUC_list)
  )
}

# ----------------------------
# UI
# ----------------------------
ui <- fluidPage(
  titlePanel("1-コンパートメントモデルでの薬物動態シミュレーター（Pharmacokinetic Simulations using 1-Compartment Model）"),
  
  tabsetPanel(
    type = "tabs",
    
    # --------------------------------------------------
    # Tab 1
    # --------------------------------------------------
    tabPanel(
      "IV Bolus Single Dose",
      sidebarLayout(
        sidebarPanel(
          textInput("drug_name1", "モデル薬物名 (Drug Name)", value = "Drug A"),
          numericInput("dose1", "投与量 (Dose)：Div (mg)", value = 500),
          numericInput("V1", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight1", "体重 (Body Weight) (kg)", value = 70),
          
          radioButtons(
            "cl_or_ke1", "計算方法を選択 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke1 == 'CLtot'",
            numericInput("CLtot1", "全身クリアランス (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke1 == 'ke'",
            numericInput("ke1", "消失速度定数 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke1 == 't1/2'",
            numericInput("t_half1", "半減期 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("time1", "シミュレーション時間 (Simulation Time) (hr)", value = 24),
          numericInput("effective_range_lower1", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper1", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot1"),
          textOutput("AUC1"),
          textOutput("results1"),
          textOutput("crossingTimes1"),
          tags$h4("薬物動態計算式："),
          withMathJax("$$Cp = \\frac{Dose}{Vd} \\times e^{-ke \\times Time}$$"),
          tags$h4("t1/2, ke, CLtotの関係式："),
          withMathJax("$$t_{1/2} = \\frac{\\ln(2)}{ke}$$"),
          withMathJax("$$ke = \\frac{\\ln(2)}{t_{1/2}}$$"),
          withMathJax("$$CLtot = ke \\times Vd$$")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 2
    # --------------------------------------------------
    tabPanel(
      "Oral Single Dose",
      sidebarLayout(
        sidebarPanel(
          textInput("drug_name2", "モデル薬物名 (Drug Name)", value = "Drug B"),
          numericInput("dose2", "投与量 (Dose)：Dpo (mg)", value = 500),
          numericInput("V2", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight2", "体重 (Body Weight) (kg)", value = 70),
          numericInput("ka2", "吸収速度定数 (ka) (1/hr)", value = 1.0),
          numericInput("F2", "バイオアベイラビリティ (F)", value = 0.8),
          
          radioButtons(
            "cl_or_ke2", "計算方法を選択 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke2 == 'CLtot'",
            numericInput("CLtot2", "全身クリアランス (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke2 == 'ke'",
            numericInput("ke2", "消失速度定数 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke2 == 't1/2'",
            numericInput("t_half2", "半減期 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("time2", "シミュレーション時間 (Simulation Time) (hr)", value = 24),
          numericInput("effective_range_lower2", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper2", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot2"),
          textOutput("AUC2"),
          textOutput("results2"),
          textOutput("crossingTimes2"),
          textOutput("CmaxTmax2"),
          tags$h4("薬物動態計算式："),
          withMathJax("$$Cp = F \\times \\frac{Dose}{Vd} \\times \\frac{ka}{ka - ke} \\times (e^{-ke \\times Time} - e^{-ka \\times Time})$$"),
          tags$h4("t1/2, ke, CLtotの関係式："),
          withMathJax("$$t_{1/2} = \\frac{\\ln(2)}{ke}$$"),
          withMathJax("$$ke = \\frac{\\ln(2)}{t_{1/2}}$$"),
          withMathJax("$$CLtot = ke \\times Vd$$")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 3
    # --------------------------------------------------
    tabPanel(
      "IV Drip Infusion",
      sidebarLayout(
        sidebarPanel(
          textInput("drug_name3", "モデル薬物名 (Drug Name)", value = "Drug C"),
          numericInput("rate3", "投与速度 (Rate)：R (mg/hr)", value = 50),
          numericInput("V3", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight3", "体重 (Body Weight) (kg)", value = 70),
          numericInput("infusion_time3", "投与時間 (Infusion Time) (hr)", value = 2),
          
          radioButtons(
            "cl_or_ke3", "計算方法を選択 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke3 == 'CLtot'",
            numericInput("CLtot3", "全身クリアランス (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke3 == 'ke'",
            numericInput("ke3", "消失速度定数 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke3 == 't1/2'",
            numericInput("t_half3", "半減期 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("time3", "シミュレーション時間 (Simulation Time) (hr)", value = 24),
          numericInput("effective_range_lower3", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper3", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot3"),
          textOutput("AUC3"),
          textOutput("results3"),
          textOutput("crossingTimes3"),
          tags$h4("薬物動態計算式："),
          tags$h5("投与中："),
          withMathJax("$$Cp = \\frac{Rate}{ke \\times Vd} \\times (1 - e^{-ke \\times Time})$$"),
          tags$h5("投与後："),
          withMathJax("$$Cp = \\frac{Rate}{ke \\times Vd} \\times (1 - e^{-ke \\times Infusion\\ Time}) \\times e^{-ke \\times (Time - Infusion\\ Time)}$$"),
          tags$h4("t1/2, ke, CLtotの関係式："),
          withMathJax("$$t_{1/2} = \\frac{\\ln(2)}{ke}$$"),
          withMathJax("$$ke = \\frac{\\ln(2)}{t_{1/2}}$$"),
          withMathJax("$$CLtot = ke \\times Vd$$")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 4
    # --------------------------------------------------
    tabPanel(
      "Oral Single Dose (CLtot=CLh+CLr)",
      sidebarLayout(
        sidebarPanel(
          tags$h4("正常時 (Normal Condition)"),
          textInput("drug_name4", "モデル薬物名 (Drug Name)", value = "Drug D"),
          numericInput("dose4", "変更前投与量 (Dose, Initial) (mg)", value = 500),
          numericInput("V4", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("height4", "身長 (Height) (cm)", value = 170),
          numericInput("weight4", "体重 (Body Weight) (kg)", value = 70),
          selectInput("sex4", "性別 (Sex)", choices = c("男性", "女性")),
          numericInput("age4", "年齢 (Age) (years)", value = 40),
          numericInput("scr4", "血清クレアチニン (Scr) (mg/dL)", value = 1),
          numericInput("ka4", "吸収速度定数 (ka) (1/hr)", value = 1.0),
          numericInput("F4", "バイオアベイラビリティ (F)", value = 0.8),
          numericInput("CLtot4", "正常時 全身クリアランス (CLtot) (L/hr)", value = 6),
          numericInput("Ae4", "尿中未変化体排泄率 (Ae)", value = 0.3),
          
          tags$hr(),
          tags$h4("病態時 (Disease Condition)"),
          numericInput("CLr_disease4", "病態時 腎クリアランス (CLr, Disease) (L/hr)", value = 2),
          numericInput("CLh_disease4", "病態時 肝クリアランス (CLh, Disease) (L/hr)", value = 2),
          numericInput("dose_disease4", "病態時 変更後投与量 (Dose, Changed) (mg)", value = 250),
          
          tags$hr(),
          numericInput("time4", "シミュレーション時間 (Simulation Time) (hr)", value = 24),
          numericInput("effective_range_lower4", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper4", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot4"),
          textOutput("BSA4"),
          textOutput("eGFR4"),
          textOutput("Ccr4"),
          textOutput("AUC4"),
          textOutput("results4"),
          textOutput("crossingTimes4"),
          tags$h4("3本の血中濃度曲線を表示しています。青：正常時、赤：病態時、緑：病態時に薬物投与量を変更した場合。"),
          tags$h4("正常時の条件をまず入力して、血中濃度曲線、CLr, CLh, AUCを確認してください。青色曲線。"),
          tags$h4("その後、病態に伴う変更内容の指示に従い、CLr、CLhがどのように変化するかを各自計算して、病態時CLr、CLhを入力してください。赤色曲線。"),
          tags$h4("正常時と病態時で各種パラメーターがどのように変化するか、確認してください。さらに、正常時と同じAUCとするためには、投与量をどのように変更すればよいか検討してください。緑色曲線。"),
          tags$h4("eGFRの計算関数（DuBois式を使用）："),
          withMathJax("$$BSA = 0.007184 \\times (Weight^{0.425}) \\times (Height^{0.725})$$"),
          withMathJax("男性: $$eGFR = 175 \\times (Scr^{-1.154}) \\times (Age^{-0.203}) \\times \\frac{BSA}{1.73}$$"),
          withMathJax("女性: $$eGFR = 175 \\times (Scr^{-1.154}) \\times (Age^{-0.203}) \\times 0.742 \\times \\frac{BSA}{1.73}$$"),
          tags$h4("Ccrの計算関数（Cockcroft-Gault式）："),
          withMathJax("男性: $$Ccr = \\frac{(140 - Age) \\times Weight}{72 \\times Scr}$$"),
          withMathJax("女性: $$Ccr = \\frac{(140 - Age) \\times Weight}{72 \\times Scr} \\times 0.85$$"),
          tags$h4("薬物動態計算式："),
          withMathJax("$$CLtot = ke \\times Vd = \\frac{F \\times Dpo}{AUC}$$"),
          withMathJax("$$CLr = CLtot \\times Ae$$"),
          withMathJax("$$CLh = CLtot - CLr$$"),
          withMathJax("$$Cp = F \\times \\frac{Dose}{Vd} \\times \\frac{ka}{ka - ke} \\times (e^{-ke \\times Time} - e^{-ka \\times Time})$$")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 5
    # --------------------------------------------------
    tabPanel(
      "IV Bolus Repeated Dose",
      sidebarLayout(
        sidebarPanel(
          textInput("drug_name5", "モデル薬物名 (Drug Name)", value = "Drug E"),
          numericInput("dose5", "投与量 (Dose)：Div (mg)", value = 500),
          numericInput("V5", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight5", "体重 (Body Weight) (kg)", value = 70),
          
          radioButtons(
            "cl_or_ke5", "計算方法を選択 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke5 == 'CLtot'",
            numericInput("CLtot5", "全身クリアランス (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke5 == 'ke'",
            numericInput("ke5", "消失速度定数 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke5 == 't1/2'",
            numericInput("t_half5", "半減期 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("tau5", "投与間隔 (τ) (hr)", value = 12),
          numericInput("n5", "投与回数 (Number of Doses)", value = 5),
          numericInput("time5", "シミュレーション時間 (Simulation Time) (hr)", value = 100),
          numericInput("effective_range_lower5", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper5", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot5"),
          textOutput("AUC5"),
          textOutput("totalAUC5"),
          textOutput("AUC24h5"),
          textOutput("results5"),
          textOutput("crossingTimes5"),
          textOutput("CssAve5"),
          tags$h4("薬物動態計算式："),
          withMathJax("$$Cp = \\sum_{i=0}^{n-1} \\left( \\frac{Dose}{Vd} \\times e^{-ke \\times (Time - i \\times \\tau)} \\times (Time \\geq i \\times \\tau) \\right)$$"),
          tags$h4("この数式は、各投与後の濃度を計算し、それを合計する形で表現しています。"),
          tags$h4("t1/2, ke, CLtotの関係式："),
          withMathJax("$$t_{1/2} = \\frac{\\ln(2)}{ke}$$"),
          withMathJax("$$ke = \\frac{\\ln(2)}{t_{1/2}}$$"),
          withMathJax("$$CLtot = ke \\times Vd$$"),
          tags$h4("＜参考＞定常状態での平均血中濃度計算式："),
          withMathJax("$$Css,ave = \\frac{Dose}{CLtot \\times \\tau}$$")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 6
    # --------------------------------------------------
    tabPanel(
      "IV Bolus Repeated, Changed Dose",
      sidebarLayout(
        sidebarPanel(
          textInput("drug_name6", "モデル薬物名 (Drug Name)", value = "Drug F"),
          numericInput("dose6_1", "初期投与量 (Initial Dose)：Div (mg)", value = 500),
          numericInput("dose6_2", "変更後の投与量 (Changed Dose)：Div (mg)", value = 250),
          numericInput("V6", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight6", "体重 (Body Weight) (kg)", value = 70),
          
          radioButtons(
            "cl_or_ke6", "計算方法を選択 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke6 == 'CLtot'",
            numericInput("CLtot6", "全身クリアランス (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke6 == 'ke'",
            numericInput("ke6", "消失速度定数 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke6 == 't1/2'",
            numericInput("t_half6", "半減期 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("tau6", "投与間隔 (τ) (hr)", value = 12),
          numericInput("n6", "総投与回数 (Total Number of Doses)", value = 5),
          numericInput("initial_dose_num6", "初期量投与回数 (Number of Initial Doses)", value = 3),
          
          numericInput("time6", "シミュレーション時間 (Simulation Time) (hr)", value = 100),
          numericInput("effective_range_lower6", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper6", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot6"),
          textOutput("AUC6"),
          textOutput("totalAUC6"),
          textOutput("AUC24h6"),
          textOutput("results6"),
          textOutput("crossingTimes6"),
          textOutput("CssAve6"),
          tags$h4("薬物動態計算式："),
          withMathJax("$$Cp = \\sum_{i=0}^{n-1} \\left( \\frac{Dose_i}{Vd} \\times e^{-ke \\times (Time - i \\times \\tau)} \\times (Time \\geq i \\times \\tau) \\right)$$"),
          tags$h4("この数式は、各投与後の濃度を計算し、それを合計する形で表現しています。"),
          tags$h4("t1/2, ke, CLtotの関係式："),
          withMathJax("$$t_{1/2} = \\frac{\\ln(2)}{ke}$$"),
          withMathJax("$$ke = \\frac{\\ln(2)}{t_{1/2}}$$"),
          withMathJax("$$CLtot = ke \\times Vd$$"),
          tags$h4("＜参考＞定常状態での平均血中濃度計算式："),
          withMathJax("$$Css,ave = \\frac{Dose}{CLtot \\times \\tau}$$")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 7
    # --------------------------------------------------
    tabPanel(
      "Oral Repeated Dose",
      sidebarLayout(
        sidebarPanel(
          textInput("drug_name7", "モデル薬物名 (Drug Name)", value = "Drug G"),
          numericInput("dose7", "投与量 (Dose)：Dpo (mg)", value = 500),
          numericInput("V7", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight7", "体重 (Body Weight) (kg)", value = 70),
          numericInput("ka7", "吸収速度定数 (ka) (1/hr)", value = 1),
          numericInput("F7", "バイオアベイラビリティ (F)", value = 0.8),
          
          radioButtons(
            "cl_or_ke7", "計算方法を選択 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke7 == 'CLtot'",
            numericInput("CLtot7", "全身クリアランス (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke7 == 'ke'",
            numericInput("ke7", "消失速度定数 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke7 == 't1/2'",
            numericInput("t_half7", "半減期 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("tau7", "投与間隔 (τ) (hr)", value = 12),
          numericInput("n7", "投与回数 (Number of Doses)", value = 5),
          numericInput("time7", "シミュレーション時間 (Simulation Time) (hr)", value = 100),
          numericInput("effective_range_lower7", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper7", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot7"),
          textOutput("AUC7"),
          textOutput("totalAUC7"),
          textOutput("AUC24h7"),
          textOutput("results7"),
          textOutput("crossingTimes7"),
          textOutput("CssAve7"),
          tags$h4("薬物動態計算式（経口繰り返し投与）："),
          withMathJax("
$$Cp = \\sum_{i=0}^{n-1} \\left(
F \\times \\frac{Dose}{Vd} \\times \\frac{ka}{ka - ke}
\\times \\left(e^{-ke \\times (Time - i \\times \\tau)} - e^{-ka \\times (Time - i \\times \\tau)}\\right)
\\times (Time \\geq i \\times \\tau)
\\right)$$
"),
          tags$h4("t1/2, ke, CLtotの関係式："),
          withMathJax("$$t_{1/2} = \\frac{\\ln(2)}{ke}$$"),
          withMathJax("$$ke = \\frac{\\ln(2)}{t_{1/2}}$$"),
          withMathJax("$$CLtot = ke \\times Vd$$"),
          tags$h4("＜参考＞定常状態での平均血中濃度計算式："),
          withMathJax("$$Css,ave = \\frac{F \\times Dose}{CLtot \\times \\tau}$$")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 8
    # --------------------------------------------------
    tabPanel(
      "Oral Repeated, 2 Comparison",
      sidebarLayout(
        sidebarPanel(
          tags$h4("条件1 (Condition 1)"),
          textInput("drug_name8_1", "モデル薬物名 1 (Drug Name)", value = "Drug H-1"),
          numericInput("dose8_1", "投与量 1 (Dose) (mg)", value = 500),
          numericInput("V8_1", "体重当たりの分布容積 1 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight8_1", "体重 1 (Body Weight) (kg)", value = 70),
          numericInput("ka8_1", "吸収速度定数 1 (ka) (1/hr)", value = 1),
          numericInput("F8_1", "バイオアベイラビリティ 1 (F)", value = 0.8),
          
          radioButtons(
            "cl_or_ke8_1", "計算方法 1 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke8_1 == 'CLtot'",
            numericInput("CLtot8_1", "全身クリアランス 1 (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke8_1 == 'ke'",
            numericInput("ke8_1", "消失速度定数 1 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke8_1 == 't1/2'",
            numericInput("t_half8_1", "半減期 1 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("tau8_1", "投与間隔 1 (τ) (hr)", value = 12),
          numericInput("n8_1", "投与回数 1 (Number of Doses)", value = 10),
          
          tags$hr(),
          
          tags$h4("条件2 (Condition 2)"),
          textInput("drug_name8_2", "モデル薬物名 2 (Drug Name)", value = "Drug H-2"),
          numericInput("dose8_2", "投与量 2 (Dose) (mg)", value = 250),
          numericInput("V8_2", "体重当たりの分布容積 2 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("weight8_2", "体重 2 (Body Weight) (kg)", value = 70),
          numericInput("ka8_2", "吸収速度定数 2 (ka) (1/hr)", value = 1),
          numericInput("F8_2", "バイオアベイラビリティ 2 (F)", value = 0.8),
          
          radioButtons(
            "cl_or_ke8_2", "計算方法 2 (Choose Calculation Method)",
            choices = list(
              "全身クリアランス (Total Clearance)：CLtot" = "CLtot",
              "消失速度定数 (Elimination Rate Constant)：ke" = "ke",
              "半減期 (Half-Life)：t1/2" = "t1/2"
            )
          ),
          conditionalPanel(
            condition = "input.cl_or_ke8_2 == 'CLtot'",
            numericInput("CLtot8_2", "全身クリアランス 2 (CLtot) (L/hr)", value = 6)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke8_2 == 'ke'",
            numericInput("ke8_2", "消失速度定数 2 (ke) (1/hr)", value = 0.1)
          ),
          conditionalPanel(
            condition = "input.cl_or_ke8_2 == 't1/2'",
            numericInput("t_half8_2", "半減期 2 (t1/2) (hr)", value = 6)
          ),
          
          numericInput("tau8_2", "投与間隔 2 (τ) (hr)", value = 24),
          numericInput("n8_2", "投与回数 2 (Number of Doses)", value = 5),
          
          tags$hr(),
          numericInput("time8", "シミュレーション時間 (Simulation Time) (hr)", value = 150),
          numericInput("effective_range_lower8", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper8", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot8"),
          textOutput("AUC8_1"),
          textOutput("totalAUC8_1"),
          textOutput("AUC24h8_1"),
          textOutput("CssAve8_1"),
          textOutput("results8_1"),
          tags$hr(),
          textOutput("AUC8_2"),
          textOutput("totalAUC8_2"),
          textOutput("AUC24h8_2"),
          textOutput("CssAve8_2"),
          textOutput("results8_2"),
          tags$hr(),
          textOutput("crossingTimes8"),
          tags$h4("薬物動態計算式（経口繰り返し投与）："),
          withMathJax("
$$Cp = \\sum_{i=0}^{n-1} \\left(
F \\times \\frac{Dose}{Vd} \\times \\frac{ka}{ka - ke}
\\times \\left(e^{-ke \\times (Time - i \\times \\tau)} - e^{-ka \\times (Time - i \\times \\tau)}\\right)
\\times (Time \\geq i \\times \\tau)
\\right)$$
"),
          tags$h4("赤：条件1（Red is for Condition 1）、青：条件2（Blue is for Condition 2）")
        )
      )
    ),
    
    # --------------------------------------------------
    # Tab 9
    # --------------------------------------------------
    tabPanel(
      "Oral Repeated (CLtot=CLh+CLr)",
      sidebarLayout(
        sidebarPanel(
          tags$h4("正常時 (Normal Condition)"),
          textInput("drug_name9", "モデル薬物名 (Drug Name)", value = "Drug I"),
          numericInput("dose9", "正常時 投与量 (Dose, Normal) (mg)", value = 500),
          numericInput("V9", "体重当たりの分布容積 (Vd/kg) (L/kg)", value = 0.6),
          numericInput("height9", "身長 (Height) (cm)", value = 170),
          numericInput("weight9", "体重 (Body Weight) (kg)", value = 70),
          selectInput("sex9", "性別 (Sex)", choices = c("男性", "女性")),
          numericInput("age9", "年齢 (Age) (years)", value = 40),
          numericInput("scr9", "血清クレアチニン (Scr) (mg/dL)", value = 1),
          numericInput("ka9", "吸収速度定数 (ka) (1/hr)", value = 1.0),
          numericInput("F9", "バイオアベイラビリティ (F)", value = 0.8),
          numericInput("CLtot9", "正常時 全身クリアランス (CLtot) (L/hr)", value = 6),
          numericInput("Ae9", "尿中未変化体排泄率 (Ae)", value = 0.3),
          numericInput("tau9", "正常時 投与間隔 (τ, Normal) (hr)", value = 24),
          numericInput("n9", "投与回数 (Number of Doses)", value = 10),
          
          
          tags$hr(),
          tags$h4("病態時 (Disease Condition)"),
          numericInput("CLr_disease9", "病態時 腎クリアランス (CLr, Disease) (L/hr)", value = 2),
          numericInput("CLh_disease9", "病態時 肝クリアランス (CLh, Disease) (L/hr)", value = 1),
          numericInput("dose_disease9", "病態時 変更後投与量 (Dose, Changed) (mg)", value = 250),
          numericInput("tau_disease9", "病態時 変更後投与間隔 (τ, Changed) (hr)", value = 24),
          numericInput("n_disease9", "病態時 変更後投与回数 (Number of Doses, Changed)", value = 10),
          
          tags$hr(),
          numericInput("time9", "シミュレーション時間 (Simulation Time) (hr)", value = 200),
          numericInput("effective_range_lower9", "有効濃度範囲 下限 (Lower Bound) (mg/L)", value = 5),
          numericInput("effective_range_upper9", "有効濃度範囲 上限 (Upper Bound) (mg/L)", value = 15)
        ),
        mainPanel(
          plotOutput("concentrationPlot9"),
          textOutput("BSA9"),
          textOutput("eGFR9"),
          textOutput("Ccr9"),
          textOutput("AUC9"),
          textOutput("AUC24h9"),
          textOutput("results9"),
          textOutput("crossingTimes9"),
          tags$h4("3本の血中濃度曲線を表示しています。青：正常時、赤：病態時、緑：病態時に投与量・投与間隔を変更した場合。"),
          tags$h4("正常時の条件をまず入力して、血中濃度曲線、CLr, CLh, AUCを確認してください。青色曲線。"),
          tags$h4("その後、病態時の CLr, CLh を入力し、病態時の投与間隔も必要に応じて変更してください。赤色曲線。"),
          tags$h4("さらに、投与量と投与間隔を変更した場合の濃度推移を確認してください。緑色曲線。"),
          tags$h4("eGFRの計算関数（DuBois式を使用）："),
          withMathJax("$$BSA = 0.007184 \\times (Weight^{0.425}) \\times (Height^{0.725})$$"),
          withMathJax("男性: $$eGFR = 175 \\times (Scr^{-1.154}) \\times (Age^{-0.203}) \\times \\frac{BSA}{1.73}$$"),
          withMathJax("女性: $$eGFR = 175 \\times (Scr^{-1.154}) \\times (Age^{-0.203}) \\times 0.742 \\times \\frac{BSA}{1.73}$$"),
          tags$h4("Ccrの計算関数（Cockcroft-Gault式）："),
          withMathJax("男性: $$Ccr = \\frac{(140 - Age) \\times Weight}{72 \\times Scr}$$"),
          withMathJax("女性: $$Ccr = \\frac{(140 - Age) \\times Weight}{72 \\times Scr} \\times 0.85$$"),
          tags$h4("薬物動態計算式（経口繰り返し投与）："),
          withMathJax("
$$Cp = \\sum_{i=0}^{n-1} \\left(
F \\times \\frac{Dose}{Vd} \\times \\frac{ka}{ka - ke}
\\times \\left(e^{-ke \\times (Time - i \\times \\tau)} - e^{-ka \\times (Time - i \\times \\tau)}\\right)
\\times (Time \\geq i \\times \\tau)
\\right)$$
"),
          tags$h4("CLtot, CLr, CLh の関係式："),
          withMathJax("$$CLr = CLtot \\times Ae$$"),
          withMathJax("$$CLh = CLtot - CLr$$"),
          withMathJax("$$ke = \\frac{CLtot}{Vd}$$")
        )
      )
    )
  )
)

# ----------------------------
# server
# ----------------------------
server <- function(input, output) {
  
  # ---------------- Tab 1 ----------------
  output$concentrationPlot1 <- renderPlot({
    drug_name <- input$drug_name1
    dose <- input$dose1
    Vd <- input$V1 * input$weight1
    Time <- seq(0, input$time1, by = 0.1)
    
    ke <- switch(
      input$cl_or_ke1,
      "CLtot" = input$CLtot1 / Vd,
      "ke" = input$ke1,
      "t1/2" = t_half_to_ke(input$t_half1)
    )
    
    sim <- simulate_iv_bolus_single(dose, Vd, ke, Time)
    Cp <- sim$Cp
    AUC <- sim$AUC
    
    crossing_times <- calc_crossing_times(Time, Cp, input$effective_range_lower1)
    crossing_text <- paste("Crossing Times (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    output$AUC1 <- renderText({
      paste("AUC for one dose:", round(AUC, 2), "mg*hr/L")
    })
    
    output$results1 <- renderText({
      paste(
        "投与量 (Dose):", input$dose1, "mg、",
        "分布容積 (Vd):", round(Vd, 2), "L、",
        "体重:", input$weight1, "kg、",
        "消失速度定数 (ke):", round(ke, 4), "1/hr、",
        "全身クリアランス (CLtot):", round(ke * Vd, 2), "L/hr、",
        "半減期 (t1/2):", round(log(2) / ke, 2), "hr"
      )
    })
    
    output$crossingTimes1 <- renderText({ crossing_text })
    
    ggplot(data.frame(Time, Cp), aes(x = Time, y = Cp)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower1, linetype = "dashed", color = "red") +
      geom_hline(yintercept = input$effective_range_upper1, linetype = "dashed", color = "red") +
      labs(
        title = paste("IV Bolus Single Dose (", drug_name, ")", sep = ""),
        x = "Time (hours)", y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
  # ---------------- Tab 2 ----------------
  output$concentrationPlot2 <- renderPlot({
    drug_name <- input$drug_name2
    dose <- input$dose2
    Vd <- input$V2 * input$weight2
    Time <- seq(0, input$time2, by = 0.1)
    
    ke <- switch(
      input$cl_or_ke2,
      "CLtot" = input$CLtot2 / Vd,
      "ke" = input$ke2,
      "t1/2" = t_half_to_ke(input$t_half2)
    )
    
    ka <- input$ka2
    F <- input$F2
    
    sim <- simulate_oral_single(dose, Vd, ke, ka, F, Time)
    Cp <- sim$Cp
    AUC <- sim$AUC
    
    crossing_times <- calc_crossing_times(Time, Cp, input$effective_range_lower2)
    crossing_text <- paste("Crossing Times (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    Cmax <- max(Cp)
    Tmax <- Time[which.max(Cp)]
    
    output$AUC2 <- renderText({
      paste("AUC for one dose:", round(AUC, 2), "mg*hr/L")
    })
    
    output$results2 <- renderText({
      paste(
        "投与量 (Dose):", dose, "mg、",
        "分布容積 (Vd):", round(Vd, 2), "L、",
        "体重:", input$weight2, "kg、",
        "吸収速度定数 (ka):", ka, "1/hr、",
        "バイオアベイラビリティ (F):", F, "、",
        "消失速度定数 (ke):", round(ke, 4), "1/hr、",
        "全身クリアランス (CLtot):", round(ke * Vd, 2), "L/hr、",
        "半減期 (t1/2):", round(log(2) / ke, 2), "hr"
      )
    })
    
    output$crossingTimes2 <- renderText({ crossing_text })
    
    output$CmaxTmax2 <- renderText({
      paste("Cmax:", round(Cmax, 2), "mg/L, Tmax:", round(Tmax, 2), "hours")
    })
    
    ggplot(data.frame(Time, Cp), aes(x = Time, y = Cp)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower2, linetype = "dashed", color = "red") +
      geom_hline(yintercept = input$effective_range_upper2, linetype = "dashed", color = "red") +
      labs(
        title = paste("Oral Single Dose (", drug_name, ")", sep = ""),
        x = "Time (hours)", y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
  # ---------------- Tab 3 ----------------
  output$concentrationPlot3 <- renderPlot({
    drug_name <- input$drug_name3
    rate <- input$rate3
    Vd <- input$V3 * input$weight3
    Time <- seq(0, input$time3, by = 0.1)
    infusion_time <- input$infusion_time3
    
    ke <- switch(
      input$cl_or_ke3,
      "CLtot" = input$CLtot3 / Vd,
      "ke" = input$ke3,
      "t1/2" = t_half_to_ke(input$t_half3)
    )
    
    sim <- simulate_iv_infusion(rate, Vd, ke, infusion_time, Time)
    Cp <- sim$Cp
    AUC <- sim$AUC
    
    crossing_times <- calc_crossing_times(Time, Cp, input$effective_range_lower3)
    crossing_text <- paste("Crossing Times (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    output$AUC3 <- renderText({
      paste("AUC:", round(AUC, 2), "mg*hr/L")
    })
    
    output$results3 <- renderText({
      paste(
        "投与速度 (Rate):", rate, "mg/hr、",
        "分布容積 (Vd):", round(Vd, 2), "L、",
        "体重:", input$weight3, "kg、",
        "投与時間:", infusion_time, "hr、",
        "消失速度定数 (ke):", round(ke, 4), "1/hr、",
        "全身クリアランス (CLtot):", round(ke * Vd, 2), "L/hr、",
        "半減期 (t1/2):", round(log(2) / ke, 2), "hr"
      )
    })
    
    output$crossingTimes3 <- renderText({ crossing_text })
    
    ggplot(data.frame(Time, Cp), aes(x = Time, y = Cp)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower3, linetype = "dashed", color = "red") +
      geom_hline(yintercept = input$effective_range_upper3, linetype = "dashed", color = "red") +
      labs(
        title = paste("IV Drip Infusion (", drug_name, ")", sep = ""),
        x = "Time (hours)", y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
  # ---------------- Tab 4 ----------------
  output$concentrationPlot4 <- renderPlot({
    drug_name <- input$drug_name4
    Time <- seq(0, input$time4, by = 0.1)
    
    dose_normal <- input$dose4
    dose_changed <- input$dose_disease4
    
    Vd <- input$V4 * input$weight4
    CLtot_normal <- input$CLtot4
    ke_normal <- CLtot_normal / Vd
    ka <- input$ka4
    F <- input$F4
    Ae <- input$Ae4
    
    CLr_normal <- CLtot_normal * Ae
    CLh_normal <- CLtot_normal - CLr_normal
    
    CLr_disease <- input$CLr_disease4
    CLh_disease <- input$CLh_disease4
    CLtot_disease <- CLr_disease + CLh_disease
    ke_disease <- CLtot_disease / Vd
    
    Cp_normal <- F * dose_normal / Vd * (ka / (ka - ke_normal)) * (exp(-ke_normal * Time) - exp(-ka * Time))
    Cp_disease <- F * dose_normal / Vd * (ka / (ka - ke_disease)) * (exp(-ke_disease * Time) - exp(-ka * Time))
    Cp_changed <- F * dose_changed / Vd * (ka / (ka - ke_disease)) * (exp(-ke_disease * Time) - exp(-ka * Time))
    
    AUC_normal <- F * dose_normal / CLtot_normal
    AUC_disease <- F * dose_normal / CLtot_disease
    AUC_changed <- F * dose_changed / CLtot_disease
    
    BSA <- 0.007184 * (input$weight4^0.425) * (input$height4^0.725)
    eGFR <- calculate_eGFR(input$scr4, input$age4, input$weight4, input$height4, input$sex4)
    Ccr <- calculate_Ccr(input$scr4, input$age4, input$weight4, input$sex4)
    
    crossing_times <- calc_crossing_times(Time, Cp_normal, input$effective_range_lower4)
    crossing_text <- paste("Crossing Times of normal curve (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    output$BSA4 <- renderText({
      paste("BSA:", round(BSA, 2), "m^2")
    })
    output$eGFR4 <- renderText({
      paste("eGFR:", round(eGFR, 2), "mL/min/1.73m^2換算ではなく体表面積補正後の推定値")
    })
    output$Ccr4 <- renderText({
      paste("Ccr:", round(Ccr, 2), "mL/min")
    })
    output$AUC4 <- renderText({
      paste(
        "AUC | 正常時:", round(AUC_normal, 2),
        "mg*hr/L, 病態時:", round(AUC_disease, 2),
        "mg*hr/L, 病態時変更後:", round(AUC_changed, 2), "mg*hr/L"
      )
    })
    output$results4 <- renderText({
      paste(
        "正常時 CLtot:", round(CLtot_normal, 2), "L/hr、",
        "正常時 ke:", round(ke_normal, 4), "1/hr、",
        "正常時 CLr:", round(CLr_normal, 2), "L/hr、",
        "正常時 CLh:", round(CLh_normal, 2), "L/hr、",
        "病態時 CLtot:", round(CLtot_disease, 2), "L/hr、",
        "病態時 ke:", round(ke_disease, 4), "1/hr、",
        "病態時 CLr:", round(CLr_disease, 2), "L/hr、",
        "病態時 CLh:", round(CLh_disease, 2), "L/hr"
      )
    })
    output$crossingTimes4 <- renderText({ crossing_text })
    
    df <- rbind(
      data.frame(Time = Time, Cp = Cp_normal, Condition = "Normal"),
      data.frame(Time = Time, Cp = Cp_disease, Condition = "Disease"),
      data.frame(Time = Time, Cp = Cp_changed, Condition = "Disease + Changed Dose")
    )
    
    ggplot(df, aes(x = Time, y = Cp, color = Condition)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower4, linetype = "dashed", color = "darkgreen") +
      geom_hline(yintercept = input$effective_range_upper4, linetype = "dashed", color = "darkgreen") +
      labs(
        title = paste("Oral Single Dose (CLtot = CLh + CLr) (", drug_name, ")", sep = ""),
        x = "Time (hours)", y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
  # ---------------- Tab 5 ----------------
  output$concentrationPlot5 <- renderPlot({
    drug_name <- input$drug_name5
    dose <- input$dose5
    Vd <- input$V5 * input$weight5
    Time <- seq(0, input$time5, by = 0.1)
    tau <- input$tau5
    n <- input$n5
    
    ke <- switch(
      input$cl_or_ke5,
      "CLtot" = input$CLtot5 / Vd,
      "ke" = input$ke5,
      "t1/2" = t_half_to_ke(input$t_half5)
    )
    CLtot <- ke * Vd
    
    sim <- simulate_iv_bolus_repeated(dose, Vd, ke, tau, n, Time)
    Cp <- sim$Cp
    AUC_list <- sim$AUC_list
    AUC_total <- sim$totalAUC
    AUC_24h <- calc_auc24(Time, Cp)
    Css_ave <- dose / (CLtot * tau)
    
    crossing_times <- calc_crossing_times(Time, Cp, input$effective_range_lower5)
    crossing_text <- paste("Crossing Times (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    output$AUC5 <- renderText({
      paste("AUC for each dose:", paste(round(AUC_list, 2), collapse = ", "), "mg*hr/L")
    })
    output$totalAUC5 <- renderText({
      paste("Total AUC:", round(AUC_total, 2), "mg*hr/L")
    })
    output$AUC24h5 <- renderText({
      paste("AUC for each 24h:", paste(round(AUC_24h, 2), collapse = ", "), "mg*hr/L")
    })
    output$results5 <- renderText({
      paste(
        "投与量:", dose, "mg、",
        "Vd:", round(Vd, 2), "L、",
        "体重:", input$weight5, "kg、",
        "ke:", round(ke, 4), "1/hr、",
        "CLtot:", round(CLtot, 2), "L/hr、",
        "t1/2:", round(log(2) / ke, 2), "hr、",
        "投与間隔:", tau, "hr、",
        "投与回数:", n, "、",
        "有効濃度範囲:", input$effective_range_lower5, "-", input$effective_range_upper5, "mg/L"
      )
    })
    output$crossingTimes5 <- renderText({ crossing_text })
    output$CssAve5 <- renderText({
      paste("定常状態での平均血中濃度 (Css, ave):", round(Css_ave, 2), "mg/L")
    })
    
    ggplot(data.frame(Time, Cp), aes(x = Time, y = Cp)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower5, linetype = "dashed", color = "red") +
      geom_hline(yintercept = input$effective_range_upper5, linetype = "dashed", color = "red") +
      labs(
        title = paste("IV Bolus Repeated Dose (", drug_name, ")", sep = ""),
        x = "Time (hours)", y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
  # ---------------- Tab 6 ----------------
  output$concentrationPlot6 <- renderPlot({
    drug_name <- input$drug_name6
    initial_dose <- input$dose6_1
    changed_dose <- input$dose6_2
    initial_dose_num <- input$initial_dose_num6
    Vd <- input$V6 * input$weight6
    Time <- seq(0, input$time6, by = 0.1)
    tau <- input$tau6
    n <- input$n6
    
    ke <- switch(
      input$cl_or_ke6,
      "CLtot" = input$CLtot6 / Vd,
      "ke" = input$ke6,
      "t1/2" = t_half_to_ke(input$t_half6)
    )
    CLtot <- ke * Vd
    
    sim <- simulate_iv_bolus_repeated_changed(
      initial_dose = initial_dose,
      changed_dose = changed_dose,
      initial_dose_num = initial_dose_num,
      Vd = Vd,
      ke = ke,
      tau = tau,
      n = n,
      Time = Time
    )
    
    Cp <- sim$Cp
    AUC_list <- sim$AUC_list
    AUC_total <- sim$totalAUC
    AUC_24h <- calc_auc24(Time, Cp)
    Css_ave <- changed_dose / (CLtot * tau)
    
    crossing_times <- calc_crossing_times(Time, Cp, input$effective_range_lower6)
    crossing_text <- paste("Crossing Times (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    output$AUC6 <- renderText({
      paste("AUC for each dose:", paste(round(AUC_list, 2), collapse = ", "), "mg*hr/L")
    })
    output$totalAUC6 <- renderText({
      paste("Total AUC:", round(AUC_total, 2), "mg*hr/L")
    })
    output$AUC24h6 <- renderText({
      paste("AUC for each 24h:", paste(round(AUC_24h, 2), collapse = ", "), "mg*hr/L")
    })
    output$results6 <- renderText({
      paste(
        "初期投与量:", initial_dose, "mg、",
        "変更後投与量:", changed_dose, "mg、",
        "総投与回数:", n, "、",
        "初期量投与回数:", initial_dose_num, "、",
        "Vd:", round(Vd, 2), "L、",
        "ke:", round(ke, 4), "1/hr、",
        "CLtot:", round(CLtot, 2), "L/hr、",
        "t1/2:", round(log(2) / ke, 2), "hr、",
        "投与間隔:", tau, "hr"
      )
    })
    output$crossingTimes6 <- renderText({ crossing_text })
    output$CssAve6 <- renderText({
      paste("変更後投与量ベースの Css, ave:", round(Css_ave, 2), "mg/L")
    })
    
    ggplot(data.frame(Time, Cp), aes(x = Time, y = Cp)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower6, linetype = "dashed", color = "red") +
      geom_hline(yintercept = input$effective_range_upper6, linetype = "dashed", color = "red") +
      labs(
        title = paste("IV Bolus Repeated, Changed Dose (", drug_name, ")", sep = ""),
        x = "Time (hours)", y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
  # ---------------- Tab 7 ----------------
  output$concentrationPlot7 <- renderPlot({
    drug_name <- input$drug_name7
    dose <- input$dose7
    Vd <- input$V7 * input$weight7
    Time <- seq(0, input$time7, by = 0.1)
    tau <- input$tau7
    n <- input$n7
    ka <- input$ka7
    F <- input$F7
    
    ke <- switch(
      input$cl_or_ke7,
      "CLtot" = input$CLtot7 / Vd,
      "ke" = input$ke7,
      "t1/2" = t_half_to_ke(input$t_half7)
    )
    CLtot <- ke * Vd
    
    sim <- simulate_oral_repeated(dose, Vd, ke, ka, F, tau, n, Time)
    Cp <- sim$Cp
    AUC_list <- sim$AUC_list
    AUC_total <- sim$totalAUC
    AUC_24h <- calc_auc24(Time, Cp)
    Css_ave <- (F * dose) / (CLtot * tau)
    
    crossing_times <- calc_crossing_times(Time, Cp, input$effective_range_lower7)
    crossing_text <- paste("Crossing Times (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    output$AUC7 <- renderText({
      paste("AUC for each dose:", paste(round(AUC_list, 2), collapse = ", "), "mg*hr/L")
    })
    output$totalAUC7 <- renderText({
      paste("Total AUC:", round(AUC_total, 2), "mg*hr/L")
    })
    output$AUC24h7 <- renderText({
      paste("AUC for each 24h:", paste(round(AUC_24h, 2), collapse = ", "), "mg*hr/L")
    })
    output$results7 <- renderText({
      paste(
        "投与量:", dose, "mg、",
        "Vd:", round(Vd, 2), "L、",
        "体重:", input$weight7, "kg、",
        "ka:", ka, "1/hr、",
        "F:", F, "、",
        "ke:", round(ke, 4), "1/hr、",
        "CLtot:", round(CLtot, 2), "L/hr、",
        "t1/2:", round(log(2) / ke, 2), "hr、",
        "投与間隔:", tau, "hr、",
        "投与回数:", n
      )
    })
    output$crossingTimes7 <- renderText({ crossing_text })
    output$CssAve7 <- renderText({
      paste("定常状態での平均血中濃度 (Css, ave):", round(Css_ave, 2), "mg/L")
    })
    
    ggplot(data.frame(Time, Cp), aes(x = Time, y = Cp)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower7, linetype = "dashed", color = "red") +
      geom_hline(yintercept = input$effective_range_upper7, linetype = "dashed", color = "red") +
      labs(
        title = paste("Oral Repeated Dose (", drug_name, ")", sep = ""),
        x = "Time (hours)", y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
  # ---------------- Tab 8 ----------------
  output$concentrationPlot8 <- renderPlot({
    Time <- seq(0, input$time8, by = 0.1)
    
    # 条件1
    drug_name1 <- input$drug_name8_1
    dose1 <- input$dose8_1
    Vd1 <- input$V8_1 * input$weight8_1
    tau1 <- input$tau8_1
    n1 <- input$n8_1
    ka1 <- input$ka8_1
    F1 <- input$F8_1
    
    ke1 <- switch(
      input$cl_or_ke8_1,
      "CLtot" = input$CLtot8_1 / Vd1,
      "ke" = input$ke8_1,
      "t1/2" = t_half_to_ke(input$t_half8_1)
    )
    CLtot1 <- ke1 * Vd1
    
    sim1 <- simulate_oral_repeated(dose1, Vd1, ke1, ka1, F1, tau1, n1, Time)
    Cp1 <- sim1$Cp
    AUC_list1 <- sim1$AUC_list
    AUC_total1 <- sim1$totalAUC
    AUC_24h_1 <- calc_auc24(Time, Cp1)
    Css_ave1 <- (F1 * dose1) / (CLtot1 * tau1)
    
    # 条件2
    drug_name2 <- input$drug_name8_2
    dose2 <- input$dose8_2
    Vd2 <- input$V8_2 * input$weight8_2
    tau2 <- input$tau8_2
    n2 <- input$n8_2
    ka2 <- input$ka8_2
    F2 <- input$F8_2
    
    ke2 <- switch(
      input$cl_or_ke8_2,
      "CLtot" = input$CLtot8_2 / Vd2,
      "ke" = input$ke8_2,
      "t1/2" = t_half_to_ke(input$t_half8_2)
    )
    CLtot2 <- ke2 * Vd2
    
    sim2 <- simulate_oral_repeated(dose2, Vd2, ke2, ka2, F2, tau2, n2, Time)
    Cp2 <- sim2$Cp
    AUC_list2 <- sim2$AUC_list
    AUC_total2 <- sim2$totalAUC
    AUC_24h_2 <- calc_auc24(Time, Cp2)
    Css_ave2 <- (F2 * dose2) / (CLtot2 * tau2)
    
    crossing_times1 <- calc_crossing_times(Time, Cp1, input$effective_range_lower8)
    crossing_times2 <- calc_crossing_times(Time, Cp2, input$effective_range_lower8)
    
    output$crossingTimes8 <- renderText({
      paste0(
        "Crossing Times (lower bound) 条件1: ",
        paste(round(crossing_times1, 2), collapse = ", "),
        " / 条件2: ",
        paste(round(crossing_times2, 2), collapse = ", ")
      )
    })
    
    output$AUC8_1 <- renderText({
      paste("条件1 AUC for each dose:", paste(round(AUC_list1, 2), collapse = ", "), "mg*hr/L")
    })
    output$totalAUC8_1 <- renderText({
      paste("条件1 Total AUC:", round(AUC_total1, 2), "mg*hr/L")
    })
    output$AUC24h8_1 <- renderText({
      paste("条件1 AUC for each 24h:", paste(round(AUC_24h_1, 2), collapse = ", "), "mg*hr/L")
    })
    output$CssAve8_1 <- renderText({
      paste("条件1 Css, ave:", round(Css_ave1, 2), "mg/L")
    })
    output$results8_1 <- renderText({
      paste(
        "条件1 | 投与量:", dose1, "mg、",
        "Vd:", round(Vd1, 2), "L、",
        "体重:", input$weight8_1, "kg、",
        "ka:", ka1, "1/hr、",
        "F:", F1, "、",
        "ke:", round(ke1, 4), "1/hr、",
        "CLtot:", round(CLtot1, 2), "L/hr、",
        "t1/2:", round(log(2) / ke1, 2), "hr、",
        "τ:", tau1, "hr、",
        "投与回数:", n1
      )
    })
    
    output$AUC8_2 <- renderText({
      paste("条件2 AUC for each dose:", paste(round(AUC_list2, 2), collapse = ", "), "mg*hr/L")
    })
    output$totalAUC8_2 <- renderText({
      paste("条件2 Total AUC:", round(AUC_total2, 2), "mg*hr/L")
    })
    output$AUC24h8_2 <- renderText({
      paste("条件2 AUC for each 24h:", paste(round(AUC_24h_2, 2), collapse = ", "), "mg*hr/L")
    })
    output$CssAve8_2 <- renderText({
      paste("条件2 Css, ave:", round(Css_ave2, 2), "mg/L")
    })
    output$results8_2 <- renderText({
      paste(
        "条件2 | 投与量:", dose2, "mg、",
        "Vd:", round(Vd2, 2), "L、",
        "体重:", input$weight8_2, "kg、",
        "ka:", ka2, "1/hr、",
        "F:", F2, "、",
        "ke:", round(ke2, 4), "1/hr、",
        "CLtot:", round(CLtot2, 2), "L/hr、",
        "t1/2:", round(log(2) / ke2, 2), "hr、",
        "τ:", tau2, "hr、",
        "投与回数:", n2
      )
    })
    
    plot_df <- rbind(
      data.frame(Time = Time, Cp = Cp1, Condition = paste0("条件1: ", drug_name1)),
      data.frame(Time = Time, Cp = Cp2, Condition = paste0("条件2: ", drug_name2))
    )
    
    ggplot(plot_df, aes(x = Time, y = Cp, color = Condition)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower8, linetype = "dashed", color = "darkgreen") +
      geom_hline(yintercept = input$effective_range_upper8, linetype = "dashed", color = "darkgreen") +
      labs(
        title = "Oral Repeated, 2 Comparison",
        x = "Time (hours)",
        y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  

  # ---------------- Tab 9 ----------------
  output$concentrationPlot9 <- renderPlot({
    drug_name <- input$drug_name9
    Time <- seq(0, input$time9, by = 0.1)
    
    # 正常時
    dose_normal <- input$dose9
    Vd <- input$V9 * input$weight9
    CLtot_normal <- input$CLtot9
    ke_normal <- CLtot_normal / Vd
    ka <- input$ka9
    F <- input$F9
    Ae <- input$Ae9
    tau_normal <- input$tau9
    
    CLr_normal <- CLtot_normal * Ae
    CLh_normal <- CLtot_normal - CLr_normal
    
    # 病態時
    CLr_disease <- input$CLr_disease9
    CLh_disease <- input$CLh_disease9
    CLtot_disease <- CLr_disease + CLh_disease
    ke_disease <- CLtot_disease / Vd
    
    dose_changed <- input$dose_disease9
    tau_disease <- input$tau_disease9
    
    # 投与回数は simulation time に基づいて自動設定
    n_normal <- floor(input$time9 / tau_normal) + 1
    n_disease <- floor(input$time9 / tau_normal) + 1
    n_changed <- floor(input$time9 / tau_disease) + 1
    
    # 濃度計算
    sim_normal <- simulate_oral_repeated(dose_normal, Vd, ke_normal, ka, F, tau_normal, n_normal, Time)
    Cp_normal <- sim_normal$Cp
    AUC_total_normal <- sim_normal$totalAUC
    AUC24_normal <- calc_auc24(Time, Cp_normal)
    
    sim_disease <- simulate_oral_repeated(dose_normal, Vd, ke_disease, ka, F, tau_normal, n_disease, Time)
    Cp_disease <- sim_disease$Cp
    AUC_total_disease <- sim_disease$totalAUC
    AUC24_disease <- calc_auc24(Time, Cp_disease)
    
    sim_changed <- simulate_oral_repeated(dose_changed, Vd, ke_disease, ka, F, tau_disease, n_changed, Time)
    Cp_changed <- sim_changed$Cp
    AUC_total_changed <- sim_changed$totalAUC
    AUC24_changed <- calc_auc24(Time, Cp_changed)
    
    # eGFR, Ccr, BSA
    BSA <- 0.007184 * (input$weight9^0.425) * (input$height9^0.725)
    eGFR <- calculate_eGFR(input$scr9, input$age9, input$weight9, input$height9, input$sex9)
    Ccr <- calculate_Ccr(input$scr9, input$age9, input$weight9, input$sex9)
    
    # crossing time
    crossing_times <- calc_crossing_times(Time, Cp_normal, input$effective_range_lower9)
    crossing_text <- paste("Crossing Times of normal curve (lower bound):", paste(round(crossing_times, 2), collapse = ", "))
    
    output$BSA9 <- renderText({
      paste("BSA:", round(BSA, 2), "m^2")
    })
    
    output$eGFR9 <- renderText({
      paste("eGFR:", round(eGFR, 2), "mL/min/1.73m^2換算ではなく体表面積補正後の推定値")
    })
    
    output$Ccr9 <- renderText({
      paste("Ccr:", round(Ccr, 2), "mL/min")
    })
    

    
    output$AUC9 <- renderText({
      paste(
        "Total AUC | 正常時:", round(AUC_total_normal, 2), "mg*hr/L, ",
        "病態時:", round(AUC_total_disease, 2), "mg*hr/L, ",
        "病態時変更後:", round(AUC_total_changed, 2), "mg*hr/L"
      )
    })
    
    output$AUC24h9 <- renderText({
      paste(
        "AUC for each 24h | 正常時:", paste(round(AUC24_normal, 2), collapse = ", "),
        " / 病態時:", paste(round(AUC24_disease, 2), collapse = ", "),
        " / 病態時変更後:", paste(round(AUC24_changed, 2), collapse = ", ")
      )
    })
    
    output$results9 <- renderText({
      paste(
        "正常時 CLtot:", round(CLtot_normal, 2), "L/hr、",
        "正常時 ke:", round(ke_normal, 4), "1/hr、",
        "正常時 CLr:", round(CLr_normal, 2), "L/hr、",
        "正常時 CLh:", round(CLh_normal, 2), "L/hr、",
        "正常時 τ:", round(tau_normal, 2), "hr || ",
        "病態時 CLtot:", round(CLtot_disease, 2), "L/hr、",
        "病態時 ke:", round(ke_disease, 4), "1/hr、",
        "病態時 CLr:", round(CLr_disease, 2), "L/hr、",
        "病態時 CLh:", round(CLh_disease, 2), "L/hr、",
        "病態時変更後 Dose:", round(dose_changed, 2), "mg、",
        "病態時変更後 τ:", round(tau_disease, 2), "hr"
      )
    })
    
    output$crossingTimes9 <- renderText({
      crossing_text
    })
    
    df <- rbind(
      data.frame(Time = Time, Cp = Cp_normal, Condition = "Normal"),
      data.frame(Time = Time, Cp = Cp_disease, Condition = "Disease"),
      data.frame(Time = Time, Cp = Cp_changed, Condition = "Disease + Changed")
    )
    
    ggplot(df, aes(x = Time, y = Cp, color = Condition)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = input$effective_range_lower9, linetype = "dashed", color = "darkgreen") +
      geom_hline(yintercept = input$effective_range_upper9, linetype = "dashed", color = "darkgreen") +
      labs(
        title = paste("Oral Repeated (CLtot=CLh+CLr) (", drug_name, ")", sep = ""),
        x = "Time (hours)",
        y = "Concentration (mg/L)"
      ) +
      theme_bw(base_size = 18)
  })
  
}

# ----------------------------
# Run App
# ----------------------------
shinyApp(ui = ui, server = server)

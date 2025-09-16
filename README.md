# HW2 - Relative Orientation (Swing-Swing Method)  
相對方位 (雙像旋轉法) [共面式解法]

本專案為 **攝影測量作業二**，目標是透過雙像旋轉法 (Swing-Swing Method)，利用共面條件進行 **解析法相對方位 (Analytical Relative Orientation, ARO)** 計算，並重建影像對的三維空間模型。

---

## 📂 專案內容

- `HW2.pdf`  
  作業報告，包含背景、公式推導、程式流程、計算成果與分析。

- `HW RO swing-swing.pdf`  
  課堂投影片與作業說明。

- `101678xy.txt`  
  像點量測座標檔 (左片/右片對應點)。

- `HW2.m`  
  MATLAB 程式碼，用於執行相對方位計算。

- `ro.f`  
  Fortran 程式碼 (老師提供版本，經過格式修改)。

- `ro.exe`  
  Fortran 編譯後的可執行檔 (Windows 下可直接執行)。

---

## 🎯 作業目的

- 了解 **共面條件 (Coplanarity condition)** 的意義與數學表達。  
- 學習 **雙像旋轉法** 的相對方位原理與計算流程。  
- 使用 **最小二乘法與迭代收斂**，解算五個相對方位元素。  
- 比較 MATLAB 與 Fortran 的程式成果，驗證運算正確性與精度。  

---

## 📸 實驗資料

- 相機率定檔 `144116.txt`：  
  - 焦距 f = 152.818 mm  
  - 主點校正值 (x0, y0)  
  - 其他內方位參數  

- 像點量測檔 `101678xy.txt`：  
  - 格式：點號、左片(x,y)、右片(x,y)、code  
  - 單位：µm (後續轉換為 mm)  
  - 共有 65 個同名點  

---

## 🛠️ 計算方法

### 1. 資料處理 (Python)
- 讀取兩張影像的測點，合併同名點。  
- 單位轉換 µm → mm。  
- 生成輸入檔 `output.txt` (格式：點號、左x、左y、右x、右y)。  

### 2. 相對方位計算 (MATLAB)
- 輸入初始值 (phie_L, kapa_L, omega_R, phie_R, kapa_R = 0°)。  
- 設定基線長 B = 40 mm，焦距 f = 152.818 mm，收斂閾值 1e-8。  
- 透過 **共面條件式線性化**，進行最小二乘與迭代。  
- 輸出：
  - 迭代結果與收斂過程  
  - 五個未知數 (相對方位角) 與其後驗中誤差  
  - 模型點座標 (U,V,W)  
  - 偽觀測體積 F、vol 與 RMS  
  - 單位權中誤差  

### 3. Fortran 程式 (ro.f / ro.exe)
- `ro.f` 為原始碼，可修改與重編譯。  
- `ro.exe` 為編譯後版本，可直接輸入資料進行計算。  
- 與 MATLAB 程式計算結果比較，驗證數值一致性。  

---

## 📊 成果與分析

- **自由度**：65 - 5 = 60  
- **五個相對方位角 (MATLAB 結果)**：  
  - Phie_L = -0.674575 ± 0.004335°  
  - Kapa_L = -2.078596 ± 0.009487°  
  - Omega_R = -0.549328 ± 0.003293°  
  - Phie_R = -0.575121 ± 0.003606°  
  - Kapa_R = -0.138761 ± 0.009536°  

- **RMS (體積均方根)**：55.999921 mm³  
- **單位權中誤差**：1.41089234 mm² (門檻 3σ ≈ 4.23 mm²)  
- **MATLAB 與 Fortran 比較**：  
  - 數值幾乎一致，差異僅在小數點末位。  
  - Phie 角在 MATLAB 為負值、Fortran 為正值，推測為座標系差異。  

---

## 📌 問題回答 (作業要求)

1. 共 65 個物點參與計算。  
2. 65 個觀測值，5 個未知數，自由度 = 60。  
3. 精度最佳元素：Omega_R (±0.003293°)。  
4. 其他角度誤差倍數：Phie_L 1.316 倍、Kapa_L 2.881 倍、Phie_R 1.095 倍、Kapa_R 2.885 倍。  
5. 模型坐標範圍：  
   - U = -23.295 ~ 75.010 mm  
   - V = -69.887 ~ 63.164 mm  
   - W = -100.541 ~ -94.657 mm  
6. 虛擬體積觀測 RMS = 55.999921 mm³  
7. 後驗單位權中誤差 = 1.41089234 mm²  
8. 意義：反映實際測量誤差水準，用以評估成果精度。  

---

## 📑 參考文獻

1. MinGW-w64 編譯器 [https://www.mingw-w64.org/](https://www.mingw-w64.org/)  
2. Fortran 官網快速入門 [https://fortran-lang.org/zh_CN/learn/quickstart/](https://fortran-lang.org/zh_CN/learn/quickstart/)  
3. 課本與課堂投影片  
4. Matlab 簡易教學 [HackMD](https://hackmd.io/@FbUJsF5qTbyirb8qlvu2Fw/Sktejk7hc)  
5. ChatGPT 說明與程式輔助  

---

## 📎 附錄

- `HW2.m`：MATLAB 相對方位計算程式  
- Python 前處理程式 (合併像點座標，生成 `output.txt`)  
- `ro.f`：Fortran 原始程式碼 (可重新編譯)  
- `ro.exe`：編譯後可執行檔 (Windows)  

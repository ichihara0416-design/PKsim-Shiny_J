# PKsim Shiny

## Overview / 概要

PKsim Shiny is an educational pharmacokinetic simulation application
based on a linear one-compartment model.

PKsim Shinyは、線形1コンパートメントモデルに基づく
薬物動態教育用シミュレーションアプリケーションです。

## Intended use / 使用目的

This application is intended for pharmacokinetics education.
It is not intended for clinical decision-making or patient care.

本アプリは薬物動態教育を目的としています。
実際の患者の診療、投与設計、臨床判断には使用しないでください。

## Requirements / 必要環境

- R
- RStudio
- shiny
- ggplot2
- markdown
- showtext
- sysfonts
- curl

## How to run / 実行方法

1. Install R and RStudio.
2. Download `PKsim_1C.R`.
3. Open the file in RStudio.
4. Click "Run App" or run `shiny::runApp()`.

1. RおよびRStudioをインストールします。
2. `PKsim_1C.R`をダウンロードします。
3. RStudioでファイルを開きます。
4. "Run App"をクリックするか、`shiny::runApp()`を実行します。

## Language / 言語

The user interface includes Japanese and English labels.
Source-code comments are mainly written in Japanese.

画面には日本語と英語の表記を含みます。
ソースコード内のコメントは主に日本語です。

## License

MIT License

## Disclaimer / 免責事項

The pharmacokinetic models and assumptions in this application are
simplified for educational and research purposes.

The authors and copyright holders provide the software without
warranty and accept no responsibility for decisions or outcomes
resulting from its use.

Do not use the application or its outputs for clinical purposes.

本アプリの薬物動態モデルおよび仮定は、教育・研究目的に
単純化されています。

本ソフトウェアは無保証で提供され、著作者および著作権者は、
その利用に基づく判断または結果について責任を負いません。

本アプリまたはその出力を臨床目的に使用しないでください。

## Installation / インストール

Required R packages are installed automatically when the application
is first run. An internet connection may be required to install
packages and download the Japanese font.

初回起動時に必要なRパッケージを自動的にインストールします。
パッケージおよび日本語フォントの取得には、
インターネット接続が必要な場合があります。
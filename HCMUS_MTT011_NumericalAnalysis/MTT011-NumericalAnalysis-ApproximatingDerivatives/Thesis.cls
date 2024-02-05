\usepackage[utf8]{inputenc}
% Set sans serif font to Arial
%\setsansfont{Arial}
% Set serifed font to Times New Roman
%\setmainfont{Times New Roman}

\LoadClassWithOptions{report}
\ProvidesClass{Thesis}
\usepackage[utf8]{vietnam}

\newcommand{\R}{\mathbb{R}}
\newcommand{\Z}{\mathhb{Z}}


\usepackage{amsmath, amsthm, amssymb,latexsym,amscd,amsfonts,enumerate}

\usepackage[top=3.5cm, bottom=3.0cm, left=3.5cm, right=2.0cm]{geometry} % căn lề theo quy chuẩn KLTN
\usepackage{indentfirst}
\parskip 2pt
\headsep=12pt
\linespread{1.45}


%\usepackage{mathptmx, mathpazo} %Fonts
%\usepackage{mathptmx} %Fonts
%\usepackage{mathptmx, mathpazo}
%\usepackage[libertine,cmintegrals,cmbraces,vvarbb]{newtxmath}


%\usepackage{cite}
\usepackage{listings}
\usepackage{soul}
\usepackage[normalem]{ulem}
\useunder{\uline}{\ul}{}
\usepackage{adjustbox}
\usepackage[final]{pdfpages}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{multirow}
\usepackage{floatrow}
\newfloatcommand{capbtabbox}{table}[][\FBwidth]
\usepackage{subcaption}

\usepackage{graphicx}
\usepackage[hidelinks, colorlinks=true, allcolors=blue]{hyperref}
\usepackage{tikz}
\usetikzlibrary{calc}

% \usepackage{algorithm}
% \usepackage{algpseudocode}
% \usepackage[ruled, lined, linesnumbered, commentsnumbered, longend]{algorithm2e}
% \usepackage[algo2e]{algorithm2e}
\usepackage{algorithm, algpseudocode, amsthm}
\floatname{algorithm}{\hfil Thuật toán}
\usepackage{mathrsfs}

\usepackage{textcomp}
\usepackage{xcolor}

\usepackage{titlesec}
\titleformat{\chapter}[display]{\bfseries \large \center}{CHƯƠNG \thechapter}{0.3em}{}[]
\titleformat{\section}{\bfseries \large}{ \thesection.}{0.3em}{}[]
\titleformat{\subsection}{\it \bfseries }{ \thesubsection.}{0.3em}{}[]
\titleformat{\subsubsection}{ \it }{ \thesubsubsection.}{0.3em}{}[]
\titlespacing{\chapter}{1em}{0.1em}{3em}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{tabularx}
\usepackage{booktabs,array,multirow}
\usepackage{longtable}
\usepackage{makecell}
\usepackage[fontsize=13pt]{scrextend}
%\changefontsizes[26pt]{13pt}
%\renewcommand{\baselinestretch}{4}  %giãn dòng
\usepackage{setspace}
%\onehalfspacing

\usepackage{etoolbox}

\usepackage[abspage,user,lastpage]{zref}

\newcommand\numRoman{
	\newpage
	\pagenumbering{roman}
	\setcounter{page}{1}
}
%Đánh số bình thường
\newcommand\pageNumber{
	\newpage
	\pagenumbering{arabic}
	\setcounter{page}{1}
}

\newcommand\contents{
	\renewcommand*\contentsname{MỤC LỤC}
	\newpage
	\phantomsection
	\linespread{1.25}
	\addcontentsline{toc}{chapter}{{\bf MỤC LỤC\rm }}
	\tableofcontents%
	\linespread{1.45}
}
\newcommand\listImages{
	\renewcommand*{\listfigurename}{\bfseries DANH MỤC CÁC HÌNH VẼ, ĐỒ THỊ}
	\newpage
	\phantomsection
	\addcontentsline{toc}{chapter}{{\bf DANH MỤC CÁC HÌNH VẼ, ĐỒ THỊ\rm }}

	{% Thêm chữ hình vẽ
		\let\oldnumberline\numberline%
		\renewcommand{\numberline}{Hình~\oldnumberline}%
		\listoffigures%
	}
}

\newcommand\listTables{
	\renewcommand*{\listtablename}{\bfseries DANH MỤC CÁC BẢNG}

	{
		\let\oldnumberline\numberline%
		\renewcommand{\numberline}{Bảng~\oldnumberline}%
		\newpage
		\phantomsection
		\addcontentsline{toc}{chapter}{{\bf DANH MỤC CÁC BẢNG \rm}}
		\listoftables%

	}
}

%Định dạng tiêu đề bảng
\usepackage[format=plain,
labelfont={bf,it},
textfont=it]{caption}

\newtheorem{defn}{Định nghĩa}
\newtheorem{theorem}{Định lý}[section]
\newtheorem{lemma}{Bổ đề}
\newtheorem{prop}{Mệnh đề}
\newtheorem{prob}{Bài toán}

\newcommand\refs{
	\chapter*{TÀI LIỆU THAM KHẢO}%
	\addcontentsline{toc}{chapter}{{\bf TÀI LIỆU THAM KHẢO}}
	\phantomsection
	\titleformat{\chapter}{\bfseries  \large}{CHƯƠNG \thechapter.}{}{}[]
	\titlespacing{\chapter}{0.1em}{0.1em}{1.1em}
	\printbibliography[keyword={Vietnam},heading=subbibliography,title={Tiếng Việt:}]
	\printbibliography[notkeyword={Vietnam},heading=subbibliography,title={Tiếng Anh:}]
}

\DeclareMathOperator*{\minimize}{minimize}
\usepackage{epstopdf}

\usepackage{listings}

% insert code python
\usepackage{tcolorbox}
\tcbuselibrary{minted,breakable,xparse,skins}

\definecolor{bg}{gray}{0.95}
\DeclareTCBListing{mintedbox}{O{}m!O{}}{%
  breakable=true,
  listing engine=minted,
  listing only,
  minted language=#2,
  minted style=default,
  minted options={%
    linenos,
    gobble=0,
    breaklines=true,
    breakafter=,,
    fontsize=\small,
    numbersep=8pt,
    #1},
  boxsep=0pt,
  left skip=0pt,
  right skip=0pt,
  left=25pt,
  right=0pt,
  top=3pt,
  bottom=3pt,
  arc=5pt,
  leftrule=0pt,
  rightrule=0pt,
  bottomrule=2pt,
  toprule=2pt,
  colback=bg,
  colframe=blue!70,
  enhanced,
  overlay={%
    \begin{tcbclipinterior}
    \fill[blue!20!white] (frame.south west) rectangle ([xshift=20pt]frame.north west);
    \end{tcbclipinterior}},
  #3}

\usepackage{epigraph}
\setlength\epigraphwidth{.8\textwidth}
\setlength\epigraphrule{0pt}

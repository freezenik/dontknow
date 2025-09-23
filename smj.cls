%%
%% This is file `smj.cls',
\def\fileversion{0.3}
\def\filename{smj}
\def\filedate{2016/03/22}
%%
%% Package `smj' to use with LaTeX2e for preparation of manuscripts
%% to be submitted to Statistical Modelling (http://smj.sagepub.com)
%% License: GPL-3
%% Copyright: (C) Arnost Komarek
%% Please report errors to smj-komarek@karlin.mff.cuni.cz
%%
%% 0.2:  support for bibTeX added
%% 0.3:  some adjustments of spacing if publish option used
%% 0.4:  some more adjustment of spacing
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{smj}[\filedate\space\fileversion\space\filename\space class by Arnost Komarek]


%%%%% Options
%%%%% ------------------------------------------------------------
\newif\if@submit
\newif\if@publish
\newif\if@shortnames

\@submittrue
\@publishfalse
\@shortnamestrue

\DeclareOption{submit}{
  \@submittrue
  \@publishfalse
}

\DeclareOption{publish}{
  \@submitfalse
  \@publishtrue
}

\DeclareOption{shortnames}{
  \@shortnamestrue
}


%%%%% Process options
%%%%% -----------------------------------------------------------
\ProcessOptions

%% Basic class and its options
\if@submit
  \LoadClass[12pt, a4paper, twoside]{article}
\fi
\if@publish
  \LoadClass[11pt, a4paper, twoside]{article}
\fi

%% Required packages
\RequirePackage{graphicx}
\RequirePackage{color}
\RequirePackage{fancyvrb}
%\RequirePackage{footmisc}
\RequirePackage{amsmath,amssymb}
\RequirePackage{calc}
%%\RequirePackage[T1]{fontenc}
\RequirePackage{url}
\IfFileExists{upquote.sty}{\RequirePackage{upquote}}{}
\definecolor{red2}{rgb}{0.9333333, 0, 0}
\definecolor{mediumblue}{rgb}{0, 0, 0.8039216}
\definecolor{magenta4}{rgb}{0.545098, 0, 0.545098}
\RequirePackage[colorlinks, linkcolor=red2, citecolor=mediumblue, urlcolor=magenta4, breaklinks=true]{hyperref}
\RequirePackage{doi}

%%% ORCID stuff
%\RequirePackage{tikz}
%\usetikzlibrary{svg.path}
%%%%\definecolor{orcidlogocol}{HTML}{A6CE39}
\definecolor{orcidlogocol}{rgb}{0.651, 0.8078, 0.2235}
%\tikzset{
%   orcidlogo/.pic={
%     \fill[orcidlogocol] svg{M256,128c0,70.7-57.3,128-128,128C57.3,256,0,198.7,0,128C0,57.3,57.3,0,128,0C198.7,0,256,57.3,256,128z};
%     \fill[white] svg{M86.3,186.2H70.9V79.1h15.4v48.4V186.2z}
%     svg{M108.9,79.1h41.6c39.6,0,57,28.3,57,53.6c0,27.5-21.5,53.6-56.8,53.6h-41.8V79.1z M124.3,172.4h24.5c34.9,0,42.9-26.5,42.9-39.7c0-21.5-13.7-39.7-43.7-39.7h-23.7V172.4z}
%     svg{M88.7,56.8c0,5.5-4.5,10.1-10.1,10.1c-5.6,0-10.1-4.6-10.1-10.1c0-5.6,4.5-10.1,10.1-10.1C84.2,46.7,88.7,51.3,88.7,56.8z};
%   }
%}
%
%\newcommand\orcidicon[1]{\href{https://orcid.org/#1}{\mbox{\scalerel*{%
%				\begin{tikzpicture}[yscale=-1,transform shape]
%				\pic{orcidlogo};
%				\end{tikzpicture}
%			}{|}}}}
%\newcommand{\ORCID}[1]{$^{\mbox{\scriptsize\orcidicon{#1}}}$}

\newcommand{\ORCID}[1]{$^{\mbox{\scriptsize, \href{https://orcid.org/#1}{\textcolor{orcidlogocol}{ORCID: #1}}}}$}
%%% end of ORCID stuff

\if@submit
  \setlength{\textwidth}{6in}
  \setlength{\textheight}{9in}
  \setlength{\oddsidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\evensidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\voffset}{(\paperheight - \textheight)/2 - 2in}
  \setlength{\parindent}{0pt}
  \setlength{\parskip}{1em}

  \renewcommand{\baselinestretch}{1.8}

  %\RequirePackage{lineno}
  %\renewcommand\thelinenumber{\scriptsize\color{blue}\arabic{linenumber}\hspace{2em}}
  %\linenumbers
\fi
\if@publish
  \setlength{\textwidth}{6in}
  \setlength{\textheight}{9in}
  \setlength{\oddsidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\evensidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\voffset}{(\paperheight - \textheight)/2 - 2in}
  \setlength{\parindent}{0pt}
  \setlength{\parskip}{1em}

  \renewcommand{\baselinestretch}{1.2}  
\fi


%% Bibliography
\if@shortnames
  \usepackage[authoryear,round]{natbib}
\else
  \usepackage[authoryear,round,longnamesfirst]{natbib}
\fi
\bibpunct{(}{)}{;}{a}{,}{,}
\bibliographystyle{smj}


%% Author, title, etc.
\newcommand{\Author}[1]{\def\@Author{#1}}
\newcommand{\AuthorRunning}[1]{\def\@AuthorRunning{#1}}
\renewcommand{\and}{{\mdseries and\ }}

\newcommand{\Affil}[1]{$^{\;\mathbf{#1}}$}

\newcommand{\Title}[1]{\def\@Title{#1}}
\newcommand{\TitleRunning}[1]{\def\@TitleRunning{#1}}

\newcommand{\CorrAddress}[1]{\def\@CorrAddress{#1}}
\newcommand{\CorrEmail}[1]{\def\@CorrEmail{#1}}
\newcommand{\CorrPhone}[1]{\def\@CorrPhone{#1}}
\newcommand{\CorrFax}[1]{\def\@CorrFax{#1}}

\newcommand{\Abstract}[1]{\def\@Abstract{#1}}
\newcommand{\Keywords}[1]{\def\@Keywords{#1}}

\newcommand{\Year}[1]{\def\@Year{#1}}
\newcommand{\Volume}[1]{\def\@Volume{#1}}
\newcommand{\Issue}[1]{\def\@Issue{#1}}
\newcommand{\DOI}[1]{\def\@DOI{#1}}

\newcommand{\SubmitDate}[1]{\def\@SubmitDate{#1}}
\newcommand{\AcceptDate}[1]{\def\@AcceptDate{#1}}

\newcounter{AffilCount}              % Counter for affiliations
\newenvironment{AffilList}{          % A  list environment for the affiliations
  \begin{list}
  {$^{\arabic{AffilCount}}$}
  {\usecounter{AffilCount}
   \setlength{\topsep}{0mm}
   \setlength{\itemsep}{0mm}
   \setlength{\parskip}{0pt}
   \setlength{\parsep}{0pt}
   \setlength{\labelsep}{5pt}
   \setlength{\labelwidth}{1cm}
   \setlength{\leftmargin}{0.3cm} 
   \rmfamily
  }
}{
  \end{list}
}                 % end environment{AfEnum}

\newcommand{\Affiliations}[1]{
  \def\@Affiliations{
    \begin{AffilList}
    #1
    \end{AffilList}
  }
}



%% Defaults
\Author{Firstname Lastname\ORCID{0000-0000-0000-0000}}
\AuthorRunning{Firstname Lastname}

\Title{Title}
\TitleRunning{Running title}

\CorrAddress{Firstname Lastname, Department, University, Street 1, ZIP City, Country}
\CorrEmail{name@university.edu}
\CorrPhone{(+1)$\;$999$\;$888$\;$777}
%%\CorrFax{(+1)$\;$999$\;$888$\;$666}

\Abstract{+++ !!! +++ Abstracts must be given. +++ !!! +++}
\Keywords{+++ !!! +++ Keywords must be given. +++ !!! +++}

\Year{2999}
\Volume{00}
\Issue{0}
\DOI{10.1177/YYYYYYYXZZZZZZZZ}

\SubmitDate{2999-12-31}
\AcceptDate{2999-12-31}


%% Title page
\renewcommand{\maketitle}{
  \thispagestyle{empty}
  \markboth{\leftline{\hspace*{2em}\@AuthorRunning}}{\rightline{\@TitleRunning\hspace{2em}}}
  \pagestyle{myheadings}

  \begin{flushleft}\sffamily\bfseries
  \LARGE\@Title

  \vspace{0.2ex}
  \Large\@Author
  %\renewcommand{\@makefnmark}{}
  %\setlength{\footnotemargin}{0pt}
  %\setlength{\skip\footins}{5em}
  %\footnote{\small Address for correspondence: \@CorrAddress. E-mail: \@CorrEmail.
  %Phone: \@CorrPhone. Fax: \@CorrFax.}
  %\setcounter{footnote}{0}
  \end{flushleft}
  
  \vspace{0.2ex}
  \@Affiliations

  \vspace{0ex}
  \rule{\textwidth}{1pt}
  \textbf{Address for correspondence:} \@CorrAddress. 

  \vspace{-1em}
  \textbf{E-mail:} \texttt{\@CorrEmail}.

  \vspace{-1em}
  \textbf{Phone:} \@CorrPhone. 

  %%\vspace{-1em}
  %%\textbf{Fax:} \@CorrFax.

  \vspace{-1em}
  \rule{\textwidth}{1pt}
  \textbf{Abstract:}\enspace \@Abstract

  \vspace{0ex}
  \rule{\textwidth}{1pt}
  \textbf{Key words:}\enspace \@Keywords

}

%% Sectioning
\renewcommand{\thesection}{\arabic{section}}
\renewcommand{\section}{%
  \@startsection{section}{1}{\z@}%
%                {-2.0ex \@plus -1ex \@minus -.2ex}%
%                {0.1ex \@plus.01ex}%
                {10pt}{1pt}      
                {\sffamily\Large\bfseries}%
}

\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}
\renewcommand{\subsection}{%
  \@startsection{subsection}{2}{\z@}%
%                {-2.0ex\@plus -1ex \@minus -.1ex}%
%                {0.1ex \@plus .01ex}%
                {1pt}{1pt}  
                {\rmfamily\large\bfseries}%
}

\renewcommand{\thesubsubsection}{\arabic{section}.\arabic{subsection}.\arabic{subsubsection}}
\renewcommand{\subsubsection}{%
  \@startsection{subsubsection}{3}{\z@}%
%                {-2.0ex\@plus -1ex \@minus -.1ex}%
%                {0.1ex \@plus .01ex}%
                {1pt}{1pt}    
                {\itshape\normalsize}%
}

%% Numbering of equations
\renewcommand{\theequation}{\thesection.\arabic{equation}}
\numberwithin{equation}{section}


%% Some commands
%%%\providecommand{\url}[1]{\texttt{\textcolor{magenta4}{#1}}}

\endinput
%%
%% End of file `smj.cls'.

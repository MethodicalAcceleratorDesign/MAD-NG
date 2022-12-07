%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%  CERN Report class (single contribution)                     %
%                                                              %
%  1.0  (03/08/2005)  MG                                       %
%    - Based on cernyrep.cls                                   %
%      Adapted to 2-column articles                            %
%                                                              %
%  Last Mod.: 16 Nov. 2006 (MG)                                %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\fileversion{1.0} 
\def\filedate{2006/11/16}
\def\docdate {2006/11/16}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cernyrep2}%
    [\filedate\space cernyrep2 class (version\space\fileversion)]

% Set flag for MathTime usage
\providecommand*\CERNMTPRO{n}
\DeclareOption{mtpro}{\renewcommand{\CERNMTPRO}{y}}

% Set hyperref flag
\providecommand*\CERNHYPER{n}
\DeclareOption{hyperref}{\renewcommand{\CERNHYPER}{y}}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions
 
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions
 
% openright option starts chapters, etc. on RH page 
\LoadClass[a4paper,10pt]{report}

\RequirePackage{ifthen,calc}

% If running in hyperref mode, do not load cite package
\ifthenelse{\equal{\CERNHYPER}{y}}
 {\RequirePackage[dvips,linktocpage,colorlinks,pagebackref]{hyperref}}
 {\RequirePackage{cite}}

%% Load MathTime support files if mtpro option given

\ifthenelse{\equal{\CERNMTPRO}{y}}
 {\usepackage[mtplusscr]{mtpro}
  \usepackage[]{mtpams}
  \usepackage[mtpcal]{mtpb}}
 {}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Running headings                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage[T1]{fontenc}
\RequirePackage{cernall}
\RequirePackage{cernunits}
\RequirePackage{multicol}

\usepackage{geometry}
\geometry{% A4sized paper
   paperwidth=210mm,
   paperheight=297mm,
   right=20mm,
   bottom=23mm,
   left=20mm,
   top=24mm,
   headsep=0mm,
   footskip=0mm,
   headheight=0mm
}

%%% Level 1
\renewcommand{\savesection}{\@startsection {section}{1}{\z@}%
                                   {-0.2cm \@plus -2ex \@minus -.2ex}%
                                   {0.1ex \@plus.2ex}%
                                   {\reset@font\normalsize\bfseries}}
%%% Level 2
\renewcommand{\savesubsection}{\@startsection{subsection}{2}{\z@}%
                                     {-0.2cm \@plus -.2ex \@minus -.2ex}%
                                     {0.1ex \@plus .2ex}%
                                     {\reset@font\normalsize\bfseries}}
%%% Level 3
\renewcommand{\savesubsubsection}{\@startsection{subsubsection}{3}{\z@}%
                                     {-0.2cm \@plus -.2ex \@minus -.2ex}%
                                     {0.1ex \@plus .2ex}%
                                     {\reset@font\normalsize\bfseries\itshape}}
%%% Level 4
\renewcommand{\paragraph}{\@startsection{paragraph}{4}{\z@}%
                                     {-0.2cm \@plus -.2ex \@minus -.2ex}%
                                     {0.1ex \@plus .2ex}%
                                     {\reset@font\normalsize\itshape}}
\let\subsubsubsection\paragraph

%% Text font is Times Roman
\renewcommand{\rmdefault}{ptm}

 
%    ****************************************
%    *               TITLE                  *
%    ****************************************
%
% This definition of \maketitle taken from article.sty, except
% it does not do a \thispagestyle{plain}.
 
\def\maketitle{\par
 \begingroup
   \def\thefootnote{\fnsymbol{footnote}}
   \def\@makefnmark{\hbox
       to 0pt{$^{\@thefnmark}$\hss}}
%   \begin{multicols}{2}[\@maketitle]
 \@maketitle
 \@thanks
 \endgroup
 \setcounter{footnote}{0}
 \let\maketitle\relax
 \let\@maketitle\relax
 \gdef\@thanks{}\gdef\@author{}\gdef\@title{}\let\thanks\relax}
 
\newlength{\titleblockheight}       % so user can change it if need be
\setlength{\titleblockheight}{4cm}
 
%%       \titlblockheight introduced at 16:52:55 on 9 Jun 1990 by JOWETT
 
\iffalse
\def\@maketitle{\vbox to \titleblockheight {\hsize\textwidth
  \linewidth\hsize \vfil \centering
  {\Large \@title \par}     % Title set in \LARGE size.
  \vskip 2ex                % Vertical space after title.
  {\large \begin{tabular}[t]{c}\@author \end{tabular}\par}
  \vfil}}
\fi
\def\@maketitle{%
\begin{center}
\Large \@title\\[14pt]
\large \begin{tabular}[t]{c}\@author \end{tabular}
\end{center}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Floats                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define caption layout %%%%%

\newcommand{\setskiptab}{%
 \setlength\abovecaptionskip{1mm}
 \setlength\belowcaptionskip{1mm}}
\newcommand{\setskipfig}{%
 \setlength\abovecaptionskip{1mm}
 \setlength\belowcaptionskip{1mm}}
\iffalse
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \begin{center}\small{#1:} #2\end{center}
  \vskip\belowcaptionskip}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\small #1: #2}%
  \ifdim \wd\@tempboxa >\hsize
   \small #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\fi
\renewenvironment{table}
               {\@float{table}[H]\setskiptab\small}
               {\end@float}
\renewenvironment{table*}
               {\@dblfloat{table}[H]\setskiptab\small}
               {\end@dblfloat}
\renewenvironment{figure}
               {\@float{figure}[H]\setskipfig}
               {\end@float}
\renewenvironment{figure*}
               {\@dblfloat{figure}[H]\setskipfig}
               {\end@dblfloat}

\ifx\@Hxfloat\@Hundef\else\expandafter\endinput\fi
\let\@Hxfloat\@xfloat
\def\@xfloat#1[{\@ifnextchar{H}{\@HHfloat{#1}[}{\@Hxfloat{#1}[}}
\def\@HHfloat#1[H]{%
\expandafter\let\csname end#1\endcsname\end@Hfloat
\vskip\intextsep\vbox\bgroup\def\@captype{#1}\parindent\z@
\ignorespaces}
\def\end@Hfloat{\egroup\vskip \intextsep}
 
% The abstract Environment
 
\renewenvironment{abstract}
{\begin{center}\large\textit{Abstract}\\[4pt]\end{center}}
{\par}
 
%    ****************************************
%    *           INITIALIZATION             *
%    ****************************************
%
% Default initializations
 
%\twocolumn
\sloppy
\parindent 5mm  \parskip4pt
\leftmargini 2em
\leftmargin\leftmargini
\leftmarginv .5em
\leftmarginvi .5em
\flushbottom

\pagestyle{empty}

\endinput


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%  CERN monograph class                                        %
%                                                              %
%  1.0  (02/08/2005)  MG                                       %
%    - Based on cernall style                                  %
%  1.01 (27/10/2005)  MG                                       %
%    - Introduce mtpro option to install Mathtime fonts        %
%  1.02 (04/11/2005)  MG                                       %
%    - openright option is nolonger default                    %
%    - make 11pt default font size                             %
%    - fancyheader style uses chapter and section titles       %
%      in smallcaps                                            %
%  1.03 (10/11/2005)  MG                                       %
%    - Make thereferences environment depend on \bibname       %
%      default: ``References''                                 %
%  1.04 (12/06/2006)  MG                                       %
%    - No addvspace in lot and lof between chapter entries     %
%  1.05 (13/06/2006)  MG                                       %
%    - Introduce BIBCHAP variable to control section level     %
%      of bibliography environment                             %
%      This is controlled by the 'bibchap' option on the       %
%      \documentclass command.                                 %
%    - write \bibname in running titles for bibliography       %
%  1.06 (16/06/2006)  MG                                       %
%    - Define ptm as \rmdefault before loading mtpro.sty       %
%  1.07 (22/06/2006)  MG                                       %
%    - Add hyperref option to load hyperref.sty                %
%  1.08 (03/07/2006)  MG                                       %
%    - Do note load cite package when hyperref is loaded       %
%  1.09 (04/07/2006)  MG                                       %
%    - Make T1 default for fontenc package                     %
%  1.10 (28/01/2007)  MG                                       %
%    - Add \myurl command to show url in hypertext mode        %
%  1.11 (09/02/2010)  MG                                       %
%    - change from mtpro.sty to mtpro2.sty for mtpro option    %
%                                                              %
%  Last Mod.: 9 February 2010 (MG)                             %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\fileversion{1.11} 
\def\filedate{2010/02/09}
\def\docdate {2010/02/09}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cernmono}%
    [\filedate\space cernmono class (version\space\fileversion)]

% Set flag for MathTime usage
\providecommand*\CERNMTPRO{n}
\DeclareOption{mtpro}{\renewcommand{\CERNMTPRO}{y}}

% Set flag for title format bibliography (default section)
\providecommand*\BIBCHAP{n}
\DeclareOption{bibchap}{\renewcommand{\BIBCHAP}{y}}

% Set hyperref flag
\providecommand*\CERNHYPER{n}
\DeclareOption{hyperref}{\renewcommand{\CERNHYPER}{y}}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions

\LoadClass[a4paper,11pt]{report}
 
\RequirePackage{ifthen,calc}

% If running in hyperref mode, do not load cite package
\ifthenelse{\equal{\CERNHYPER}{y}}
 {\RequirePackage[dvips,linktocpage,colorlinks,pagebackref]{hyperref}
  \providecommand\myurl[1]{\newline[\url{#1}]\par}}
 {\RequirePackage{cite}}

%% Load MathTime support files if mtpro option given

% Declare Times Roman as default text font (should be before
% calling mtpro to get math operator font set correctly)
%    (mtpro.sty defines math operator font in function of \rmdefault)
\ifthenelse{\equal{\CERNMTPRO}{y}}
 {\usepackage[mtplusscr,mtpcal]{mtpro2}}
 {}
\iffalse
\ifthenelse{\equal{\CERNMTPRO}{y}}
 {\usepackage[mtplusscr]{mtpro}
  \usepackage[]{mtpams}
  \usepackage[mtpcal]{mtpb}}
 {}
\fi

%% Common definitions

\RequirePackage{cernall}
\RequirePackage{cernunits}
\RequirePackage[T1]{fontenc}

%%% Text font : Times Roman
\renewcommand{\rmdefault}{ptm}
 
%%% Customization of visible document elements

\renewcommand\thechapter      {\@arabic\c@chapter}
\renewcommand\thesection      {\thechapter.\@arabic\c@section}

\renewcommand\baselinestretch{1.}

\renewcommand\appendix{\par
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{figure}{0}%
  \setcounter{table}{0}%
  \setcounter{equation}{0}%
  \gdef\thesection{Appendix~\@Alph\c@section:}%
  \gdef\thesubsection{\@Alph\c@section.\@arabic\c@subsection}
  \gdef\theequation{\@Alph\c@section.\@arabic\c@equation}
  \@addtoreset {equation}{section}%
  \gdef\thetable{\@Alph\c@section.\@arabic\c@table}%
  \@addtoreset {table}{section}%
  \gdef\thefigure{\@Alph\c@section.\@arabic\c@figure}%
  \@addtoreset {figure}{section}%
}

% 
%  ********************
%  *   Part command   *
%  ********************
%

\renewcommand\part{%
  \if@openright
    \cleardoublepage
  \else
    \clearpage
  \fi
  \thispagestyle{empty}%
  \if@twocolumn
    \onecolumn
    \@tempswatrue
  \else
    \@tempswafalse
  \fi
  \null\vfil
  \secdef\@part\@spart}
                                                                                
% 
%  ************************************
%  *   thebibliography environement   *
%  ************************************
%

\setlength\bibindent{1.5em}
\renewenvironment{thebibliography}[1]
     {% First check for level of title (chapter or section)
      \ifthenelse{\equal{\BIBCHAP}{y}}
      {\chapter*{\bibname}%
       \addcontentsline{toc}{chapter}{\bibname}}%
      {\section*{\bibname}%
       \addcontentsline{toc}{section}{\bibname}}%
      \@mkboth{\bibname}{\bibname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \Nospacing
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\renewcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty

%  ****************************************
%  *            PAGE LAYOUT               *
%  ****************************************
%
 
\textheight 24cm       % Height of text (including footnotes and figures,
                       % excluding running head and foot).
\textwidth 16cm        % Width of text line.
\topmargin -28pt       % offset of top header from voffset (1 inch)
\headheight 14pt       % offset of top header from voffset (1 inch)
%\topmargin -48pt       % offset of top header from voffset (1 inch)

%%% pagewidth           = 210 mm = 596 pt
%%% hoffset :           =   1 in =  72 pt
%%% textwidth           = 160 mm = 454 pt

% SIDE MARGINS:  %%%% A4 : pagewidth=160mm 
\if@twoside   % Values for two-sided printing:
      %%% margins : inner     =  30 mm =  86 pt
      %%%           outer     =  20 mm =  56 pt
      %%% \odd/evensidemargin are with respect to hoffset -->
   \oddsidemargin   14pt    %   Left margin on odd-numbered pages.
   \evensidemargin -14pt    %   Left margin on even-numbered pages.
\else         % Values for one-sided printing:
   \oddsidemargin    0pt    %   Note that \oddsidemargin = \evensidemargin
   \evensidemargin   0pt
\fi

% PARAGRAPH INDENTATION
\setlength{\parindent}{1cm}
\setlength{\parskip}{3pt}

%% Chapter headings (when used should be like title)
\def\@chapter[#1]#2{\ifnum \c@secnumdepth >\m@ne
                         \refstepcounter{chapter}%
                         \typeout{\@chapapp\space\thechapter.}%
                         \addcontentsline{toc}{chapter}%
                                   {\protect\numberline{\thechapter}#1}%
                    \else
                      \addcontentsline{toc}{chapter}{#1}%
                    \fi
                    \chaptermark{#1}%
%MG 12062006         \addtocontents{lof}{\protect\addvspace{10\p@}}%
%MG 12062006         \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makechapterhead{#2}]%
                    \else
                      \@makechapterhead{#2}%
                      \@afterheading
                    \fi}
\def\@makechapterhead#1{%
  {\parindent \z@ \raggedright \normalfont
    \ifnum \c@secnumdepth >\m@ne
        \LARGE\bfseries \@chapapp\space \thechapter
        \par\nobreak
        \vskip 10\p@
    \fi
    \interlinepenalty\@M
    \LARGE \bfseries \boldmath #1\unboldmath\par\nobreak
    \vskip 10\p@
  }}
\def\@makeschapterhead#1{%
  {\parindent \z@ \raggedright \normalfont
    \interlinepenalty\@M
    \LARGE  \bfseries \boldmath #1\unboldmath\par\nobreak
    \vskip 10\p@
  }}

%%% make footnote robust (if it appears in argument of, e.g., \section)
%\let\oldfootnote\footnote
%\def\footnote{\protect\oldfootnote}

\iffalse
%%% Level 1
\renewcommand{\savesection}{\@startsection {section}{1}{\z@}%
                                   {-0.1pt \@plus -.1ex \@minus -.1ex}%
                                   {0.1pt \@plus.1ex}%
                                   {\reset@font\normalsize\bfseries}}
\renewcommand{\sectionstar}[1]{%\typeout{ **** #1 *****}%
                             \savesection*{\large\upshape\boldmath #1\unboldmath}}
\renewcommand{\sectionnormal}[2][!?!]{\ifthenelse
                 {\equal{#1}{!?!}}{%\typeout{ @@@@@@ section #2 @@@@@@@}%
                                   \savesection[{\upshape #2}]%
                              {\large\upshape\boldmath #2\unboldmath}}
                                  {%\typeout{ @@@@@@ section [#1]#2 @@@@@@@}%
                                   \savesection[{\upshape #1}]%
                              {\large\upshape\boldmath #2\unboldmath}}}
\def\section{\@ifstar{\sectionstar}{\sectionnormal}}

%%% Level 2
\renewcommand{\savesubsection}{\@startsection{subsection}{2}{\z@}%
                                     {-0.1pt \@plus -.1ex \@minus -.1ex}%
                                     {0.1pt \@plus .1ex}%
                                     {\reset@font\normalsize\bfseries}}
\renewcommand{\subsectionstar}[1]{%\typeout{ **** #1 *****}%
                             \savesubsection*{\boldmath #1\unboldmath}}
\renewcommand{\subsectionnormal}[2][!?!]{\ifthenelse
                 {\equal{#1}{!?!}}{%\typeout{ @@@@@@ section #2 @@@@@@@}%
                                   \savesubsection[#2]%
                              {\boldmath #2\unboldmath}}
                                  {%\typeout{ @@@@@@ section [#1]#2 @@@@@@@}%
                                   \savesubsection[#1]%
                              {\boldmath #2\unboldmath}}}
\def\subsection{\@ifstar{\subsectionstar}{\subsectionnormal}}

%%% Level 3
\renewcommand{\savesubsubsection}{\@startsection{subsubsection}{3}{\z@}%
                                     {-0.1pt \@plus -.1ex \@minus -.1ex}%
                                     {0.1pt \@plus .1ex}%
                                     {\reset@font\normalsize\bfseries\itshape}}
\renewcommand{\subsubsectionstar}[1]{%\typeout{ **** #1 *****}%
                             \savesubsubsection*{\boldmath #1\unboldmath}}
\renewcommand{\subsubsectionnormal}[2][!?!]{\ifthenelse
                 {\equal{#1}{!?!}}{%\typeout{ @@@@@@ section #2 @@@@@@@}%
                                   \savesubsubsection[#2]%
                              {\boldmath #2\unboldmath}}
                                  {%\typeout{ @@@@@@ section [#1]#2 @@@@@@@}%
                                   \savesubsubsection[#1]%
                              {\boldmath #2\unboldmath}}}
\def\subsubsection{\@ifstar{\subsubsectionstar}{\subsubsectionnormal}}

%%% Level 4
\renewcommand{\paragraph}{\@startsection{paragraph}{4}{\z@}%
                                     {-0.1pt \@plus -.1ex \@minus -.1ex}%
                                     {0.1pt \@plus .1ex}%
                                     {\reset@font\normalsize\itshape}}
\let\subsubsubsection\paragraph

%%% Level 5
% Fifth level generates run-in title
\renewcommand{\subparagraph}{\@startsection{subparagraph}{5}{\parindent}%
                                       {0.1cm \@plus.1ex \@minus .1ex}%
                                       {-1em}%
                                      {\reset@font\normalsize\bfseries}}

% The value of the counter secnumdepth gives the depth of the
% highest-level sectioning command that is to produce section numbers.
\setcounter{secnumdepth}{4}
        
\fi

\usepackage{fancyhdr}
\usepackage{truncate}
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\sectionmark}[1]{\markright{\thesection~#1}}
\renewcommand{\chaptermark}[1]{\markboth{\thechapter~#1}{}}
\renewcommand\headrulewidth{0pt}
\fancyfoot[C]{\thepage}
\fancyhead[LO]{\truncate{.94\headwidth}{\textsc{\rightmark}}}
%\fancyhead[LO]{\hspace*{3mm}\truncate{.94\headwidth}{\textsc{\rightmark}}}
%\fancyhead[lo]{\textsc{\rightmark}}
%\fancyhead[RE]{\textsc{\leftmark}\hspace*{3mm}}
\fancyhead[RE]{\textsc{\leftmark}}
\def\cleardoublepage{\clearpage\if@twoside \ifodd\c@page\else
    \thispagestyle{empty}\hbox{}\newpage
    \if@twocolumn\thispagestyle{empty}\hbox{}\newpage\fi\fi\fi}
\endinput

  \def\urlbrk{\nobreak\hskip .5em plus 1fill \penalty\z@ \hskip-.17em
  plus -1fill minus .11 em \vadjust{}\nobreak \hskip 0pt plus 1fill}
  \providecommand\myurl[2][]
  {\ifhmode\unskip\unskip\fi\textcolor{blue}{%
  \mathchardef\UrlBigBreakPenalty 7000
  \mathchardef\UrlBreakPenalty 9000
  \urlbrk#1\url{#2}}}
  }

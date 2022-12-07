%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%  CERN Report class (single contribution)                     %
%                                                              %
%  2.0  (03/08/2005)  MG                                       %
%    - Divide cernyrep.cls into cernall.sty and CERN Report    %
%      class cernyrep.cls (new)                                %
%  2.01 (27/10/2005)  MG                                       %
%    - Introduce mtpro option to install Mathtime fonts        %
%  2.02 (03/07/2006)  MG                                       %
%    - Add hyperref option to load hyperref.sty                %
%    - Do note load cite package when hyperref is loaded       %
%  2.03 (04/07/2006)  MG                                       %
%    - Make T1 default for fontenc package                     %
%  2.04 (28/01/2007)  MG                                       %
%    - Add \myurl command to show url in hypertext mode        %
%  2.05 (05/06/2007)  MG                                       %
%    - Add package ifpdf to test whether we run latex          %
%      or pdflatex, since hyperref has to be called with the   %
%      correct driver (dvips or pdftex)                        %
%  2.06 (13/05/2009)  MG                                       %
%    - Introduce option biblatex to not load cite package      %
%      and take care of some other interference between        %
%      biblatex and cernall.sty (Protect commands in           %
%      \AtBeginDocument when BIBLATEX flag set)                %
%                                                              %
%  Last Mod.: 13 May 2009 (MG)                                 %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\fileversion{2.06} 
\def\filedate{2009/05/13}
\def\docdate {2009/05/13}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cernyrep}%
    [\filedate\space cernyrep class (version\space\fileversion)]

% Set flag for MathTime usage
\providecommand*\CERNMTPRO{n}
\DeclareOption{mtpro}{\renewcommand{\CERNMTPRO}{y}}

% Set hyperref flag
\providecommand*\CERNHYPER{n}
\DeclareOption{hyperref}{\renewcommand{\CERNHYPER}{y}}

% Set biblatex flag
\providecommand*\BIBLATEX{n}
\DeclareOption{biblatex}{\renewcommand{\BIBLATEX}{y}}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions
 
% openright option starts chapters, etc. on RH page 
\LoadClass[a4paper,11pt,openright]{report}

\RequirePackage{ifthen,calc}
\RequirePackage{ifpdf}

% If running in hyperref mode, do not load cite package
\ifthenelse{\equal{\CERNHYPER}{y}}
 {\ifpdf
   \RequirePackage[pdftex,linktocpage,colorlinks,pagebackref]{hyperref}
  \else
   \RequirePackage[dvips,linktocpage,colorlinks,pagebackref]{hyperref}
  \fi
  \providecommand\myurl[1]{\newline[\url{#1}]\par}
 }
 {\ifthenelse{\equal{\BIBLATEX}{y}}%%% All these gymnastics seem to be necessary for biblatex
   {%
    \RequirePackage[UKenglish]{babel}
    \RequirePackage[babel]{csquotes}
    \RequirePackage[bibstyle=numeric-cern,citestyle=numeric-comp,block=ragged]{biblatex}
    \AtBeginDocument{\renewcommand{\figurename}{Fig.}%
                     \renewcommand{\bibname}{References}%
                     \defbibheading{bibliography}{\relax}%
                     \renewenvironment*{thebibliography}
                       {\section*{\bibname}%
                        \list
                         {\printfield[labelnumberwidth]{labelnumber}}
                         {\setlength{\labelwidth}{\labelnumberwidth}%
                          \setlength{\leftmargin}{\labelwidth}%
                          \setlength{\labelsep}{\biblabelsep}%
                          \addtolength{\leftmargin}{\labelsep}%
                          \setlength{\itemsep}{\bibitemsep}%
                          \setlength{\parsep}{\bibparsep}%
                          \renewcommand*{\makelabel}[1]{\hss##1}}%
                         \BNospacing
                         \sloppy
                         \clubpenalty4000
                         \@clubpenalty \clubpenalty
                         \widowpenalty4000%
                         \sfcode`\.\@m}
                       {\endlist}
                    }   
   }
   {\RequirePackage{cite}}%
}
%% Load MathTime support files if mtpro option given

\ifthenelse{\equal{\CERNMTPRO}{y}}
 {\usepackage[mtplusscr]{mtpro2}}
 {}
\iffalse
\ifthenelse{\equal{\CERNMTPRO}{y}}
 {\usepackage[mtplusscr]{mtpro}
  \usepackage[]{mtpams}
  \usepackage[mtpcal]{mtpb}}
 {}
\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Running headings                                                % (Removed!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \usepackage{fancyhdr}
% % \usepackage{truncate}
% \pagestyle{fancy}
% \fancyhf{}
% \renewcommand{\chaptermark}[1]{\markboth{{\scshape #1}}{}}
% \renewcommand\headrulewidth{0pt}
% \cfoot{\thepage}
% \fancyhead[LO]{\truncate{.95\headwidth}{\leftmark}}
% % The `E' option of the fancyhead package should be used in conjunction with
% % the `twoside' parameter of the \documentclass declaration. Otherwise we get
% % a warning: Package Fancyhdr Warning: \fancyhead's `E' option without twoside
% % option is useless.
% %\fancyhead[RE]{\textsc{Chapter\ \thechapter}}
% \fancyhead[R]{\textsc{Chapter\ \thechapter}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Chapter headings (when used should be like title)               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@makechapterhead#1{%
  \vspace*{1\p@}%
  {\parindent \z@ \raggedright \normalfont
    \ifnum \c@secnumdepth >\m@ne
        \Large\bfseries \@chapapp\space \thechapter
        \par\nobreak
        \vskip 15\p@
    \fi
    \interlinepenalty\@M
    \Large \bfseries \boldmath #1\unboldmath\par\nobreak
    \vskip 15\p@
  }}
\def\@makeschapterhead#1{%
  {\parindent \z@ \raggedright \normalfont
    \interlinepenalty\@M
    \LARGE  \bfseries \boldmath #1\unboldmath\par\nobreak
    \vskip 10\p@
  }}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Table of contents layout                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand*\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    \setlength\@tempdima{2em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}


%% Common definitions

\RequirePackage[T1]{fontenc}
\RequirePackage{cernall}
\RequirePackage{cernunits}

%%% Text font : Times Roman
\renewcommand{\rmdefault}{ptm}

% Display dates using the UK English format
\usepackage[UKenglish]{babel}
\usepackage[UKenglish,cleanlook]{isodate}

\endinput
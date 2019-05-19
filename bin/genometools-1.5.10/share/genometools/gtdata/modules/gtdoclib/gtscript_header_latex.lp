\chapter{gtscript Reference}
This chapter describes the \texttt{gtscript} API. \texttt{gtscript} is basically Lua (an embeddable scripting language) plus parts of the \emph{GenomeTools} C libraries exported to Lua. Because the \emph{GenomeTools} binary \texttt{gt} contains an embedded Lua interpreter, \texttt{gtscript} files can be executed with the \texttt{gt} binary. The parts of the \emph{GenomeTools} C libraries exported to Lua (the gtscript API) are described here, for a documentation of Lua itself and its APIs, please refer to the Lua reference manual.

See the index for an alphabetical list of all available interfaces.

\section*{Notes}

\begin{itemize}
    \item You have to add \texttt{require 'gtlua'} to your script in order to load the parts of the  \texttt{gtscript} API which are implemented in  \texttt{gtscript} itself.
    \item By default, all functions of the  \texttt{gtscript} API are contained in the \texttt{gt} table. That is, you have to prepend all functions calls with ``\texttt{gt}.''. For example, write \texttt{gt.show()} to call the \texttt{show()} function documented below. You can use \texttt{gt.export()} or \texttt{gt.re()} to export the \texttt{gt} table to the global environment, which makes the prepending unnecessary (but clutters your global environment).
    \item If a function is documented as ``returns an array'', this means that the function returns a table where only consecutive integer keys from 1 on are used and which can be traversed with the \texttt{ipairs()} function.
\end{itemize}

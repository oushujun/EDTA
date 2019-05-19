\section*{Table of Modules}
\label{sec:table_of_modules}
\begin{itemize}
<%for _, module in ipairs(modules) do%>
  \item Module \texttt{<%=module%>}\dotfill page \pageref{<%=module%>}
<%end%>
\end{itemize}

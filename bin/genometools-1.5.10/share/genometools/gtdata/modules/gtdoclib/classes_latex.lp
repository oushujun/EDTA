\section*{Table of Classes}
\label{sec:table_of_classes}
\begin{itemize}
<%for _, class in ipairs(classes) do%>
  \item Class \texttt{<%=class%>}\dotfill page \pageref{sec:class_<%=class%>}
<%end%>
\end{itemize}

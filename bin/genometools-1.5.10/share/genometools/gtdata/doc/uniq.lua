print ([[

A depth-first traversal of a feature node graph starts at the top-level feature
node (or pseudo-node) and explores as far along each branch as possible before
backtracking. Let's assume that the feature nodes are stored in a list in the
order of their traversal (called the ``feature node list'').

Two feature node graphs are considered to be repeated if their feature node list
(from the depth-first traversal) have the same length and each feature node pair
(from both lists at the same position) is ``similar''.

Two feature nodes are ``similar'', if they have the same sequence ID, feature
type, range, strand, and phase.

For such a repeated feature node graph the one with the higher score (of the
top-level feature) is kept. If only one of the feature node graphs has a defined
score, this one is kept.]])

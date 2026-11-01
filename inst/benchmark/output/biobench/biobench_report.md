# Secondary BioBench regression benchmark

Nine complete-case regression datasets; 30 replicates at 30% MCAR; four methods; one thread.
Runtime, pooled NRMSE/PFC and downstream known-beta coefficient absolute bias are available in
`biobench_all_results.csv`. The CSV records the seed and error text for every attempted cell.
Failed cells are excluded from dataset means and rankings, and remain visible in the raw file.

Cells represented: 1080; attempted: 1080; failed attempts: 30; not attempted: 0.
Results-selected five (metric wins, mean rank, BioBench order): anorexia, birthwt, epilepsy, theoph, phosphate.

See `biobench_selection_ranking.csv` for all nine datasets and the generated order, and
`biobench_selected_five.csv` for the selected subset.


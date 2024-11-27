<!--

Many thanks for contributing to EDTA!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

-->

## PR checklist

For Nextflow implementation,

- [ ] `conda` and `container` directives are included for each process
- [ ] Docker container + singularity container (optional) are included for each process
- [ ] Flow `meta.id` with each data channel
- [ ] Use nf-core resource labels such as `process_high`
- [ ] Used nf-core module
- [ ] Use `versions.yml` or versions topic
- [ ] No process in the `main.nf`. We can have a process in a sub-workflow file

# v0.1.0 release checklist

Do not create the release until the publication branch has been reviewed and merged into `main`.

1. Confirm the GitHub Actions unit-test workflow passes on `main`.
2. Perform the documented fresh-clone installation and synthetic example run.
3. Confirm `dota-seq-analyzer --version` reports `0.1.0` and canonical outputs report schema `3.0.0`.
4. Run the full unit-test suite and build fresh source/wheel artifacts from the tagged source.
5. Create and push an annotated tag without rewriting history:

   ```bash
   git tag -a v0.1.0 -m "DoTA-Seq Analyzer v0.1.0"
   git push origin v0.1.0
   ```

6. Create the GitHub Release from `v0.1.0`, include release notes from `CHANGELOG.md`, and attach freshly built artifacts if distributed.
7. After the repository is public and the release is final, enable the GitHub repository in Zenodo and archive the GitHub Release.
8. Add the Zenodo DOI to the release metadata and manuscript citation only after Zenodo assigns it.

The software release version (`0.1.0`) and canonical result schema version (`3.0.0`) are intentionally independent.

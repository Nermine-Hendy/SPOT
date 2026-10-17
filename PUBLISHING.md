# Publish SPOT on GitHub

## Source layout

Use one repository named `SPOT`, with `matlab/` and `python/` directories and the main README at its root. Unzip `SPOT_GitHub_Ready.zip` before publishing. The extracted source files belong in the repository; ZIPs can be attached separately to releases.

## Publish using GitHub Desktop

1. Install [GitHub Desktop](https://desktop.github.com/) and sign in to the intended account.
2. Choose **File → New repository**, name it `SPOT`, and choose a local parent folder. Leave automatic README, Git ignore and license generation disabled because the package includes them. Create the repository.
3. Open the new repository folder in the file manager. Copy the **contents** of the extracted `SPOT` folder into it, including `.gitignore` and `.gitattributes`. The result should be `SPOT/matlab/`, `SPOT/python/` and `SPOT/README.md`, without an extra nested `SPOT` folder.
4. Review the changes in GitHub Desktop. Enter `Prepare SPOT MATLAB and Python source release` as the summary, then choose **Commit to main**.
5. Choose **Publish repository**. Confirm the owner and name. Suggested description: `Optical imaging emulation for space domain awareness, with MATLAB and Python implementations.`
6. Clear **Keep this code private** when ready to make the source public, then select **Publish Repository**.
7. Use **View on GitHub** to check the main README, both implementation READMEs and folders. Add topics such as `space-domain-awareness`, `optical-imaging`, `simulation`, `matlab`, `python` and `synthetic-data` in the About settings.

The shareable URL is `https://github.com/YOUR-USERNAME/SPOT`, using the actual owner. No repository URL exists until publication succeeds.

GitHub browser uploads have a 25 MiB per-file limit. Each stellar CSV is about 41.1 MiB, so use Desktop or Git. All source-package files are below GitHub's 100 MiB regular-Git limit.

## Create a tagged release

1. On GitHub, select **Releases → Draft a new release**.
2. Choose an unused tag. `v7.2.0` is an example aligned with the MATLAB software identifier. Python identifies itself as `7.2.0-py.1`; the numerical schema is `7.0.0`. Keep repository tags and implementation identifiers distinct when reporting results.
3. Target `main`; use the title `SPOT — MATLAB and Python`.
4. Describe the included examples, dataset generation, GUI, diagnostics and test suite. State that observation-conditioned validation requires separate inputs, and that the native SatSim notebook also requires SatSim/TensorFlow.
5. Optionally attach the source ZIP. GitHub also creates source archives from the tag. Attach observation data only after confirming redistribution terms.
6. Review the tag, assets and description, then select **Publish release**.

Share the project URL for general access, `/releases/latest` for the latest release, or `/releases/tag/v7.2.0` for that specific example tag. In a paper, cite the actual tag used for the results. Once published, add the actual `repository-code` URL to `CITATION.cff`; add version/date and a publication DOI only when confirmed.

## Optional validation data

`SPOT_Validation_Inputs.zip` holds the FITS/WCS/MAT research inputs. Two FITS files exceed 100 MiB each; keep them outside regular Git commits. For local validation, copy its `Validation/` contents into `python/Inputs/Validation/`. Ignore rules protect that directory from accidental data commits.

Redistribution terms are not established by the project materials. Confirm attribution and permission before publishing the archive as a release asset or data deposit. Core simulations do not need it.

## Release checks

- Run `SPOT.check` in MATLAB with its required toolboxes.
- Run `python -m pytest -q` from `python/` in the configured environment.
- Execute the examples and validation workflows supporting the published results, retaining configuration and provenance.
- Keep notebook outputs, virtual environments, caches and generated acquisitions out of source commits.

## Official references

- [Create a repository with GitHub Desktop](https://docs.github.com/en/desktop/overview/creating-your-first-repository-using-github-desktop)
- [Publish an existing project](https://docs.github.com/en/desktop/adding-and-cloning-repositories/adding-an-existing-project-to-github-using-github-desktop)
- [File-size limits](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
- [Manage releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
- [Link to releases](https://docs.github.com/en/repositories/releasing-projects-on-github/linking-to-releases)

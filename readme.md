### Instruction for directories and files.

#### gitignore
Set in the `file.gitignore`.
Not update the data file and raw data analysis output.


### Commands for Git
steps for sync:
```bash
git status
git add uppmax_test.txt
git status
git commit -m "Add a test file from UPPMAX"
git push origin main
```

---

delete file:
```bash
git rm file
git status
git commit -m"remove file"
git push orgin main

#if the file is already deleted by rm:
git status
git add -A
git commit -m "Remove test.txt"
git push origin main
```
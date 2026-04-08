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
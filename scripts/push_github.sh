#!/bin/bash
# 将本仓库推送到 GitHub（首次使用需先登录）
set -euo pipefail
cd "$(dirname "$0")/.."

GH="${PWD}/bin/gh"
if [[ ! -x "$GH" ]]; then
  echo "请先下载 gh: https://cli.github.com/"
  exit 1
fi

if ! "$GH" auth status &>/dev/null; then
  echo "请先在终端运行: $GH auth login"
  echo "手机端可选用 GitHub 一次性验证码登录"
  exit 1
fi

REPO_NAME="${1:-meep_optimizer}"
if ! "$GH" repo view "$REPO_NAME" &>/dev/null; then
  "$GH" repo create "$REPO_NAME" --public --source=. --remote=origin --push
else
  git push -u origin main
fi

echo "完成: $("$GH" repo view --json url -q .url)"

#!/usr/bin/env bash
# setup-android.sh — Install Android SDK, NDK, and Rust cross-compilation tools
# for building harp-hymnal and harp-trainer as Android APKs.
#
# Target: aarch64-linux-android (ARM64 tablet, Android 16 / API 36)
# Host:   x86_64-unknown-linux-gnu (Ubuntu)
#
# Usage:  bash setup-android.sh
# Then:   source ~/.bashrc   (or restart terminal)

set -euo pipefail

ANDROID_HOME="${ANDROID_HOME:-$HOME/Android/Sdk}"
CMDLINE_TOOLS_VERSION="11076708"  # latest as of 2026
NDK_VERSION="27.2.12479018"
BUILD_TOOLS_VERSION="35.0.0"
PLATFORM_VERSION="android-35"  # API 35 (Android 16 runs 35+ apps fine)

echo "=== Android Cross-Compilation Setup ==="
echo "ANDROID_HOME: $ANDROID_HOME"
echo ""

# ── 1. Java (needed by Android SDK tools) ──
if ! command -v java &>/dev/null; then
    echo ">>> Installing OpenJDK 17..."
    sudo apt-get update -qq
    sudo apt-get install -y -qq openjdk-17-jdk-headless
fi
echo "Java: $(java -version 2>&1 | head -1)"

# ── 2. Android command-line tools ──
mkdir -p "$ANDROID_HOME/cmdline-tools"
if [ ! -d "$ANDROID_HOME/cmdline-tools/latest" ]; then
    echo ">>> Downloading Android command-line tools..."
    TMPZIP=$(mktemp /tmp/android-cmdline-XXXXX.zip)
    curl -fsSL -o "$TMPZIP" \
        "https://dl.google.com/android/repository/commandlinetools-linux-${CMDLINE_TOOLS_VERSION}_latest.zip"
    unzip -q "$TMPZIP" -d "$ANDROID_HOME/cmdline-tools"
    mv "$ANDROID_HOME/cmdline-tools/cmdline-tools" "$ANDROID_HOME/cmdline-tools/latest"
    rm "$TMPZIP"
fi

SDKMANAGER="$ANDROID_HOME/cmdline-tools/latest/bin/sdkmanager"

# ── 3. Accept licenses & install SDK components ──
echo ">>> Installing SDK platform, build-tools, and NDK..."
yes | "$SDKMANAGER" --licenses >/dev/null 2>&1 || true
"$SDKMANAGER" --install \
    "platform-tools" \
    "platforms;${PLATFORM_VERSION}" \
    "build-tools;${BUILD_TOOLS_VERSION}" \
    "ndk;${NDK_VERSION}"

NDK_HOME="$ANDROID_HOME/ndk/$NDK_VERSION"

# ── 4. Shell environment ──
EXPORTS=$(cat <<'ENVEOF'

# Android SDK/NDK (added by harp/setup-android.sh)
export ANDROID_HOME="$HOME/Android/Sdk"
export ANDROID_NDK_HOME="$ANDROID_HOME/ndk/27.2.12479018"
export PATH="$ANDROID_HOME/cmdline-tools/latest/bin:$ANDROID_HOME/platform-tools:$PATH"
ENVEOF
)

RCFILE="$HOME/.bashrc"
if ! grep -q "setup-android.sh" "$RCFILE" 2>/dev/null; then
    echo "$EXPORTS" >> "$RCFILE"
    echo ">>> Added ANDROID_HOME/NDK_HOME to $RCFILE"
fi

# Export for this session
export ANDROID_HOME
export ANDROID_NDK_HOME="$NDK_HOME"
export PATH="$ANDROID_HOME/cmdline-tools/latest/bin:$ANDROID_HOME/platform-tools:$PATH"

# ── 5. Rust Android target + cargo-ndk ──
echo ">>> Adding Rust target aarch64-linux-android..."
rustup target add aarch64-linux-android

if ! cargo install --list | grep -q "cargo-ndk"; then
    echo ">>> Installing cargo-ndk..."
    cargo install cargo-ndk
fi

# ── 6. Verify ──
echo ""
echo "=== Setup Complete ==="
echo "ANDROID_HOME:     $ANDROID_HOME"
echo "ANDROID_NDK_HOME: $ANDROID_NDK_HOME"
echo "NDK:              $NDK_VERSION"
echo "Rust target:      aarch64-linux-android"
echo "cargo-ndk:        $(cargo ndk --version 2>/dev/null || echo 'installed')"
echo ""
echo "Next: run ./build-android.sh to build APKs"

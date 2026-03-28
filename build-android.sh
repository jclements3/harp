#!/usr/bin/env bash
# build-android.sh — Cross-compile harp-hymnal and harp-trainer as Android APKs
#
# Prerequisites: run setup-android.sh first
# Target: aarch64-linux-android (ARM64 Android 16 tablet)
#
# Usage:
#   bash build-android.sh              # build both apps (debug)
#   bash build-android.sh --release    # build both apps (release, smaller + faster)
#   bash build-android.sh hymnal       # build only harp-hymnal
#   bash build-android.sh trainer      # build only harp-trainer

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# ── Config ──
ANDROID_HOME="${ANDROID_HOME:-$HOME/Android/Sdk}"
NDK_VERSION="27.2.12479018"
ANDROID_NDK_HOME="${ANDROID_NDK_HOME:-$ANDROID_HOME/ndk/$NDK_VERSION}"
BUILD_TOOLS="$ANDROID_HOME/build-tools/35.0.0"
PLATFORM_JAR="$ANDROID_HOME/platforms/android-35/android.jar"
TARGET="aarch64-linux-android"
MIN_SDK=28   # Android 9+ (for NativeActivity + Oboe audio)
TARGET_SDK=35
OUT_DIR="$SCRIPT_DIR/target/android-apk"

PROFILE="debug"
CARGO_PROFILE_FLAG=""
BUILD_HYMNAL=true
BUILD_TRAINER=true
BUILD_SAMPLER=true
BUILD_PLAYER=true

for arg in "$@"; do
    case "$arg" in
        --release) PROFILE="release"; CARGO_PROFILE_FLAG="--release" ;;
        hymnal)    BUILD_TRAINER=false; BUILD_SAMPLER=false; BUILD_PLAYER=false ;;
        trainer)   BUILD_HYMNAL=false; BUILD_SAMPLER=false; BUILD_PLAYER=false ;;
        sampler)   BUILD_HYMNAL=false; BUILD_TRAINER=false; BUILD_PLAYER=false ;;
        player)    BUILD_HYMNAL=false; BUILD_TRAINER=false; BUILD_SAMPLER=false ;;
    esac
done

# ── Verify tools ──
check_tool() {
    if [ ! -f "$1" ] && ! command -v "$1" &>/dev/null; then
        echo "ERROR: $1 not found. Run setup-android.sh first."
        exit 1
    fi
}

check_tool "$BUILD_TOOLS/aapt2"
check_tool "$BUILD_TOOLS/apksigner"
check_tool "$BUILD_TOOLS/zipalign"
check_tool "$PLATFORM_JAR"
check_tool cargo-ndk

if ! rustup target list --installed | grep -q "$TARGET"; then
    echo "ERROR: Rust target $TARGET not installed. Run: rustup target add $TARGET"
    exit 1
fi

mkdir -p "$OUT_DIR"

# ── Build function ──
build_apk() {
    local CRATE_DIR="$1"
    local CRATE_NAME="$2"
    local LIB_NAME="$3"
    local APP_LABEL="$4"
    local PACKAGE="$5"

    echo ""
    echo "========================================="
    echo "  Building $APP_LABEL ($PROFILE)"
    echo "========================================="

    # Step 1: Compile native library with cargo-ndk
    echo ">>> Compiling Rust -> $TARGET..."
    cd "$SCRIPT_DIR/$CRATE_DIR"
    cargo ndk \
        --target arm64-v8a \
        --platform $MIN_SDK \
        $CARGO_PROFILE_FLAG \
        -- build --lib --features android

    # Find the .so
    local SO_PATH="$SCRIPT_DIR/$CRATE_DIR/target/$TARGET/$PROFILE/lib${LIB_NAME}.so"
    if [ ! -f "$SO_PATH" ]; then
        echo "ERROR: Expected .so not found at $SO_PATH"
        echo "Looking for alternatives..."
        find "$SCRIPT_DIR/$CRATE_DIR/target/$TARGET/$PROFILE/" -name "*.so" 2>/dev/null
        exit 1
    fi
    echo ">>> Built: $SO_PATH ($(du -h "$SO_PATH" | cut -f1))"

    # Step 2: Create APK structure
    local WORK="$OUT_DIR/${CRATE_NAME}-work"
    rm -rf "$WORK"
    mkdir -p "$WORK/lib/arm64-v8a"

    # Copy .so as libmain.so (required by NativeActivity)
    cp "$SO_PATH" "$WORK/lib/arm64-v8a/libmain.so"

    # Step 3: Generate AndroidManifest.xml
    cat > "$WORK/AndroidManifest.xml" <<MANIFEST
<?xml version="1.0" encoding="utf-8"?>
<manifest xmlns:android="http://schemas.android.com/apk/res/android"
    package="$PACKAGE"
    android:versionCode="1"
    android:versionName="1.0">

    <uses-sdk android:minSdkVersion="$MIN_SDK" android:targetSdkVersion="$TARGET_SDK" />

    <!-- Microphone for pitch detection (harp-trainer) -->
    <uses-permission android:name="android.permission.RECORD_AUDIO" />

    <!-- Keep screen on during practice -->
    <uses-permission android:name="android.permission.WAKE_LOCK" />

    <!-- OpenGL ES 3.0 for egui/glow -->
    <uses-feature android:glEsVersion="0x00030000" android:required="true" />

    <application
        android:label="$APP_LABEL"
        android:icon="@mipmap/ic_launcher"
        android:roundIcon="@mipmap/ic_launcher_round"
        android:hasCode="false"
        android:debuggable="$([ "$PROFILE" = "debug" ] && echo true || echo false)">

        <activity
            android:name="android.app.NativeActivity"
            android:label="$APP_LABEL"
            android:configChanges="orientation|screenSize|screenLayout|keyboardHidden"
            android:exported="true"
            android:screenOrientation="unspecified">

            <meta-data
                android:name="android.app.lib_name"
                android:value="main" />

            <intent-filter>
                <action android:name="android.intent.action.MAIN" />
                <category android:name="android.intent.category.LAUNCHER" />
            </intent-filter>
        </activity>
    </application>
</manifest>
MANIFEST

    # Step 4: Build the APK
    echo ">>> Packaging APK..."

    local UNSIGNED_APK="$WORK/${CRATE_NAME}-unsigned.apk"
    local ALIGNED_APK="$WORK/${CRATE_NAME}-aligned.apk"
    local FINAL_APK="$OUT_DIR/${CRATE_NAME}.apk"

    # Compile icon resources if the crate has a res/ directory
    local RES_DIR="$SCRIPT_DIR/$CRATE_DIR/res"
    local LINK_RES=""
    if [ -d "$RES_DIR" ]; then
        local RES_COMPILED="$WORK/compiled_res"
        mkdir -p "$RES_COMPILED"
        "$BUILD_TOOLS/aapt2" compile --dir "$RES_DIR" -o "$RES_COMPILED/res.zip"
        LINK_RES="-R $RES_COMPILED/res.zip"
    fi

    "$BUILD_TOOLS/aapt2" link \
        -o "$UNSIGNED_APK" \
        --manifest "$WORK/AndroidManifest.xml" \
        -I "$PLATFORM_JAR" \
        $LINK_RES \
        --min-sdk-version $MIN_SDK \
        --target-sdk-version $TARGET_SDK

    # Add native libraries to the APK (it's a zip)
    cd "$WORK"
    zip -r "$UNSIGNED_APK" lib/

    # Align
    "$BUILD_TOOLS/zipalign" -f 4 "$UNSIGNED_APK" "$ALIGNED_APK"

    # Sign with debug key (auto-generated if missing)
    local KEYSTORE="$HOME/.android/debug.keystore"
    if [ ! -f "$KEYSTORE" ]; then
        mkdir -p "$HOME/.android"
        keytool -genkeypair \
            -keystore "$KEYSTORE" \
            -storepass android \
            -alias androiddebugkey \
            -keypass android \
            -keyalg RSA \
            -keysize 2048 \
            -validity 10000 \
            -dname "CN=Debug,O=Android,C=US"
    fi

    "$BUILD_TOOLS/apksigner" sign \
        --ks "$KEYSTORE" \
        --ks-pass pass:android \
        --key-pass pass:android \
        --ks-key-alias androiddebugkey \
        --out "$FINAL_APK" \
        "$ALIGNED_APK"

    # Cleanup work dir
    rm -rf "$WORK"

    echo ">>> APK ready: $FINAL_APK ($(du -h "$FINAL_APK" | cut -f1))"
}

# ── Build apps ──
if $BUILD_HYMNAL; then
    build_apk "harp-hymnal" "harp-hymnal" "harp_hymnal" "Harp Hymnal" "com.harp.hymnal"
fi

if $BUILD_TRAINER; then
    build_apk "harp-trainer" "harp-trainer" "harp_trainer" "Harp Trainer" "com.harp.trainer"
fi

if $BUILD_SAMPLER; then
    build_apk "harp-sampler" "harp-sampler" "harp_sampler" "Harp Sampler" "com.harp.sampler"
fi

if $BUILD_PLAYER; then
    build_apk "harp-player" "harp-player" "harp_player" "Harp Player" "com.harp.player"
fi

echo ""
echo "========================================="
echo "  Done! APKs in: $OUT_DIR/"
echo "========================================="
ls -lh "$OUT_DIR"/*.apk 2>/dev/null

echo ""
echo "To install on your tablet (USB connected):"
echo "  adb install $OUT_DIR/harp-hymnal.apk"
echo "  adb install $OUT_DIR/harp-trainer.apk"
echo "  adb install $OUT_DIR/harp-sampler.apk"
echo "  adb install $OUT_DIR/harp-player.apk"
echo ""
echo "Or copy the .apk files to the tablet and open them."

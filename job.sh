#!/bin/bash
#PBS -q AL
#PBS -l select=1:ncpus=1:mem=16gb
#PBS -m abe
#PBS -M fumiya@rcnp.osaka-u.ac.jp
#PBS -N garfield_sim
echo "start batch job"
# --- 環境設定 ---
#source $HOME/local/root_v6.28.06/root/bin/thisroot.sh
source /np1a/phanes/opt/root/6.32.06/bin/thisroot.sh
# For Garfield++ 2025/07/30
export CMAKE_PREFIX_PATH=$HOME/Geant4/geant4/install:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=$HOME/Geant4/qt-everywhere-src-5.15.2/qt-install/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=$HOME/Geant4/qt-everywhere-src-5.15.2/qt-install:$CMAKE_PREFIX_PATH
export PKG_CONFIG_PATH=$HOME/Geant4/qt-everywhere-src-5.15.2/qt-install/lib/pkgconfig:$PKG_CONFIG_PATH  
export CXXFLAGS="-march=native"
export CFLAGS="-march=native"
# For Geant4 2025/07/30
export PATH=$HOME/local/gcc-10.3.0/bin:$PATH
export LD_LIBRARY_PATH=$HOME/local/gcc-10.3.0/lib64:$LD_LIBRARY_PATH
source $HOME/Geant4/geant4/install/bin/geant4.sh
# For Qt6 2025/07/30 ninja
export PATH=$HOME/local:$PATH
export PATH=$HOME/Geant4/qt-everywhere-src-6.6.2/qt-install/bin:$PATH
export LD_LIBRARY_PATH=$HOME/Geant4/qt-everywhere-src-6.6.2/qt-install/lib:$LD_LIBRARY_PATH
export QT_PLUGIN_PATH=$HOME/Geant4/qt-everywhere-src-6.6.2/qt-install/plugins
export QML2_IMPORT_PATH=$HOME/Geant4/qt-everywhere-src-6.6.2/qt-install/qml
# For GCC 10.3.0
export LD_LIBRARY_PATH=$HOME/local/gcc-10.3.0/lib64:$LD_LIBRARY_PATH
# For Garfield++
export PATH=$HOME/garfieldpp-master/install:$PATH
export GARFIELD_HOME=~/garfieldpp-master/install
export LD_LIBRARY_PATH=$GARFIELD_HOME/lib64:$LD_LIBRARY_PATH 
# --- 作業ディレクトリへ移動 ---
cd $PBS_O_WORKDIR
umask 0002
# --- 日付の取得（変数に格納） ---
NOW=$(date +%Y%m%d_%H%M%S)
LOGDIR="logs/${NOW}"
# --- ログディレクトリ作成 ---
mkdir -p "$LOGDIR"
echo "check1: log dir = $LOGDIR"
# --- ログ出力 ---
echo "[INFO] Job started at $(date)" > "$LOGDIR/job.log"
echo "[INFO] Working dir: $PBS_O_WORKDIR" >> "$LOGDIR/job.log"
echo "[INFO] Node list: $(cat $PBS_NODEFILE)" >> "$LOGDIR/job.log"
# --- Garfield++ 実行 ---
./fieldview_wirecath  >> "$LOGDIR/run_sim.out" 2>> "$LOGDIR/run_sim.err"
RUN_STATUS=$?
# --- 実行完了ログ ---
if [ $RUN_STATUS -eq 0 ]; then
  echo "[INFO] Job completed successfully at $(date)" >> "$LOGDIR/job.log"
else
  echo "[ERROR] run_sim failed with exit code $RUN_STATUS at $(date)" >> "$LOGDIR/job.log"
fi
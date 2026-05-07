// visualize_tangent.C
// 円の接線探索アルゴリズム（反復法）の可視化
//
// [見どころ]
// 反復(Iteration)が進むごとに、
// 1. 直線（赤線）の角度が変わる
// 2. それに合わせて、接点（赤丸）が円周上を滑るように移動する
// 様子を確認すること。

#include "TCanvas.h"
#include "TH2D.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TF1.h"
#include "TSystem.h"
#include "TLatex.h"
#include <vector>
#include <cmath>
#include <iostream>

// 構造体: 1つのヒットデータ
struct Hit {
    double x, y; // ワイヤー位置
    double r;    // ドリフト半径
};

void visualize_tangent() {
    // 1. 仮想データの準備 (右上がりのトラック)
    std::vector<Hit> hits = {
        {-0.5, -0.2, 0.15},
        {-0.2,  0.0, 0.12},
        { 0.1,  0.2, 0.14},
        { 0.4,  0.5, 0.13}
    };

    // 2. キャンバスの準備
    TCanvas *c1 = new TCanvas("c1", "Tangent Fit Visualization", 800, 600);
    // 座標系 (-1.0 ~ 1.0)
    TH2D *frame = new TH2D("frame", "Iterative Tangent Fit;X;Y", 100, -0.8, 0.8, 100, -0.6, 1.0);
    frame->SetStats(0);

    // 3. 初期推定 (わざとズラした直線からスタート)
    // 本当は右上がりだが、最初は「水平 (a=0)」から始めてみる
    double current_a = 0.0; 
    double current_b = 0.1; // 少し上にズラす

    // --- 反復ループ (アニメーション) ---
    const int max_iter = 8;
    for (int iter = 0; iter < max_iter; ++iter) {
        
        // 描画のリセット
        c1->Clear();
        frame->Draw();

        // ---------------------------------------------
        // A. ドリフト円を描画 (固定)
        // ---------------------------------------------
        for (const auto& h : hits) {
            TEllipse *el = new TEllipse(h.x, h.y, h.r, h.r);
            el->SetFillStyle(0);
            el->SetLineColor(kGray+2);
            el->SetLineStyle(2); // 点線
            el->Draw();

            TMarker *m = new TMarker(h.x, h.y, 3); // ワイヤー中心
            m->Draw();
        }

        // ---------------------------------------------
        // B. 接点の計算 (ここが数学の心臓部)
        // ---------------------------------------------
        std::vector<double> x_tan, y_tan;
        
        // 法線ベクトルの成分 (上向き)
        double denom = std::sqrt(1.0 + current_a * current_a);
        double nx = -current_a / denom;
        double ny =  1.0       / denom;

        for (const auto& h : hits) {
            // 直線とワイヤーの距離 (符号付き)
            double dist = (current_a * h.x - h.y + current_b) / denom;
            
            // 符号の決定 (直線が上なら+1, 下なら-1)
            double sign = (dist > 0) ? 1.0 : -1.0;

            // 接点の計算（直線の角度によって位置が変わる）
            double tx = h.x + sign * h.r * nx;
            double ty = h.y + sign * h.r * ny;

            x_tan.push_back(tx);
            y_tan.push_back(ty);

            // 接点を描画 (赤丸)
            TMarker *mt = new TMarker(tx, ty, 20);
            mt->SetMarkerColor(kRed);
            mt->SetMarkerSize(1.5);
            mt->Draw();
        }

        // ---------------------------------------------
        // C. 直線の描画 (現在のパラメータ)
        // ---------------------------------------------
        TF1 *fLine = new TF1("fLine", "[0]*x + [1]", -1, 1);
        fLine->SetParameters(current_a, current_b);
        fLine->SetLineColor(kRed);
        fLine->SetLineWidth(2);
        fLine->Draw("same");

        // テキスト表示
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.15, 0.85, Form("Iteration: %d", iter));
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.15, 0.80, Form("Line: y = %.3fx + %.3f", current_a, current_b));
        latex.DrawLatex(0.15, 0.75, "Red Points = Calculated Tangent Points");

        // 画面更新
        c1->Update();
        
        // 1秒待機 (動きを見やすくするため)
        gSystem->Sleep(1000); 

        // ---------------------------------------------
        // D. 更新 (次のループのために直線をフィットし直す)
        // ---------------------------------------------
        TGraph gTemp(x_tan.size(), x_tan.data(), y_tan.data());
        // 今ある赤丸に対して直線を引く
        auto res = gTemp.Fit("pol1", "Q S N"); 
        if (res.Get()) {
            current_b = res->Parameter(0);
            current_a = res->Parameter(1);
        }
        
        // --> ここで x_tan, y_tan は破棄される (箱を空にする)
        //     次のループでは、新しい角度に基づいて再計算される
    }
}
#include "DxLib.h"
#include <math.h>
#include <algorithm>
#define PI 3.14159265358979
#define N 20
int pp(int a) {
	if (a < 0)a*=-1;
	if (a % 2 == 0)return 1;
	if (a % 2 == 1)return -1;
}
int Abs(int a){
		if (a < 0)a *= -1;
		return a;
}
double *srt[1000];
double Score = 0;
int counter = 0;

//消失点だよ
struct vpoint {
	double x = 320, y = 240;
};
typedef struct vpoint vpt;
vpt vp;
//xyzと拡大率から成る三次元の頂点だよ
struct object {
	double z = 0, x = 0, y = 0;
	double v = 0;
	int hand[10];
	int pat = 0;//障害物のパターン　0…壊せるブロック　1…壊せぬブロック　2…無敵化パネル　3…ボムHP回復パネル 自機の状態　0…常態　1…ダメージ　2…加速無敵　3…ボム　4…死んだ　
	int active = 0;//障害物が出現中かどうか 0…OFF 1…ON
	double zm = 0; //拡大率　※直接は使わない

};
typedef struct object obj;
obj pl;

//与えられたx,y,zと消失点から一点透視図法による擬似3D座標を計算するよ
obj prs(obj *org, vpt *vp) {
	obj aft;
	aft.zm = 0.00001*0.5*(org->z*org->z - org->z) + 0.0003;
	aft.x = (org->x - vp->x)*aft.zm + vp->x;
	aft.y = (org->y - vp->y)*aft.zm + vp->y;
	aft.z = org->z;
	return aft;
}
obj prs2(obj org, vpt *vp) {
	obj aft;
	aft.zm = 0.00001*0.5*(org.z*org.z - org.z) + 0.0003;
	aft.x = (org.x - vp->x)*aft.zm + vp->x;
	aft.y = (org.y - vp->y)*aft.zm + vp->y;
	aft.z = org.z;
	return aft;
}

//DrawLineを使って四角形を描くよ
void DrawRect(
	double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4, unsigned int Cr, int t) {
	//左上　右上　右下　左下
	DrawLine(x1, y1, x2, y2, Cr, t);
	DrawLine(x2, y2, x3, y3, Cr, t);
	DrawLine(x3, y3, x4, y4, Cr, t);
	DrawLine(x4, y4, x1, y1, Cr, t);
}

//与えられた座標から相対座標を計算して自機を描くよ
void DrawPlayer(obj *org, double *ang) {
	obj top = *org, rl = *org, rs = *org, rr = *org, lw1 = *org, lw2 = *org, rw1 = *org, rw2 = *org;
	obj bst[4];
	bst[0] = *org;	bst[1] = *org;	bst[2] = *org;	bst[3] = *org;
	int i;
	double L = 25;
	unsigned int Cr = GetColor(255, 155, 55);

	Cr = GetColor(155, 55, 55);
	top.z -= 2.5*L;
	rs.z += 2.5 * L;
	rl.z += 2.5 * L;
	rr.z += 2.5 * L;
	rw1.z += 3 * L;
	lw1.z += 3 * L;
	rw2.z += 1.5 * L;
	lw2.z += 1.5 * L;
	rs.x += L * cos(PI / 180 * (*ang + 180));
	rs.y += L * sin(PI / 180 * (*ang + 180));
	rl.x += L * cos(PI / 180 * (*ang + 90));
	rl.y += L * sin(PI / 180 * (*ang + 90));
	rr.x += L * cos(PI / 180 * (*ang + 270));
	rr.y += L * sin(PI / 180 * (*ang + 270));
	lw1.x += 3 * L * cos(PI / 180 * (*ang + 90));
	lw1.y += 3 * L * sin(PI / 180 * (*ang + 90));
	rw1.x += 3 * L * cos(PI / 180 * (*ang - 90));
	rw1.y += 3 * L * sin(PI / 180 * (*ang - 90));
	lw2.x = rl.x; lw2.y = rl.y;
	rw2.x = rr.x; rw2.y = rr.y;

	if (org->pat == 2) {
		SetDrawBlendMode(DX_BLENDMODE_ADD, 22 * (1 + cos(PI / 180 * org->zm * 30)) + 205);
	}
	DrawModiGraph(
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&rl, &vp).x, prs(&rl, &vp).y,
		prs(&rr, &vp).x, prs(&rr, &vp).y, org->hand[0], FALSE);
	DrawModiGraph(
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&rl, &vp).x, prs(&rl, &vp).y,
		prs(&rs, &vp).x, prs(&rs, &vp).y, org->hand[0], FALSE);
	DrawModiGraph(
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&rs, &vp).x, prs(&rs, &vp).y,
		prs(&rr, &vp).x, prs(&rr, &vp).y, org->hand[0], FALSE);

	DrawModiGraph(
		prs(&lw2, &vp).x, prs(&lw2, &vp).y,
		prs(&lw2, &vp).x, prs(&lw2, &vp).y,
		prs(&lw1, &vp).x, prs(&lw1, &vp).y,
		prs(&rl, &vp).x, prs(&rl, &vp).y, org->hand[1], FALSE);

	DrawModiGraph(
		prs(&rw2, &vp).x, prs(&rw2, &vp).y,
		prs(&rw2, &vp).x, prs(&rw2, &vp).y,
		prs(&rw1, &vp).x, prs(&rw1, &vp).y,
		prs(&rr, &vp).x, prs(&rr, &vp).y, org->hand[1], FALSE);


	DrawTriangle(
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&rl, &vp).x, prs(&rl, &vp).y,
		prs(&rr, &vp).x, prs(&rr, &vp).y, Cr, 0);
	DrawTriangle(
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&rl, &vp).x, prs(&rl, &vp).y,
		prs(&rs, &vp).x, prs(&rs, &vp).y, Cr, 0);
	DrawTriangle(
		prs(&top, &vp).x, prs(&top, &vp).y,
		prs(&rs, &vp).x, prs(&rs, &vp).y,
		prs(&rr, &vp).x, prs(&rr, &vp).y, Cr, 0);
	Cr = GetColor(255, 55, 55);
	DrawTriangle(
		prs(&rw2, &vp).x, prs(&rw2, &vp).y,
		prs(&rw1, &vp).x, prs(&rw1, &vp).y,
		prs(&rr, &vp).x, prs(&rr, &vp).y, Cr, 0);
	DrawTriangle(
		prs(&lw2, &vp).x, prs(&lw2, &vp).y,
		prs(&lw1, &vp).x, prs(&lw1, &vp).y,
		prs(&rl, &vp).x, prs(&rl, &vp).y, Cr, 0);
	if (org->pat == 2) {
		for (i = 0; i<4; i++) {
			bst[i].z += 2.5 * L;
			bst[i].x += (5 * (1 + sin(PI / 180 * org->zm * 30)) + 25) * cos(PI / 180 * (*ang + 90 * i));
			bst[i].y += (5 * (1 + sin(PI / 180 * org->zm * 30)) + 25) * sin(PI / 180 * (*ang + 90 * i));

		}
		SetDrawBlendMode(DX_BLENDMODE_ADD, 22 * (1 + cos(PI / 180 * org->zm * 30)) + 205);
		DrawModiGraph(
			prs(&bst[0], &vp).x, prs(&bst[0], &vp).y,
			prs(&bst[1], &vp).x, prs(&bst[1], &vp).y,
			prs(&bst[2], &vp).x, prs(&bst[2], &vp).y,
			prs(&bst[3], &vp).x, prs(&bst[3], &vp).y, org->hand[3], TRUE);
		SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 0);
	}
}

//与えられた座標から相対座標を計算して障害物を描くよ
void DrawObject0(obj *org, double *objang) {
	obj p1 = *org, p2 = *org, p3 = *org, p4 = *org, pr = *org, pf = *org;
	int L = 100;
	double del = 0;
	if (org->z > 480) {
		del = 3 * (org->z - 480);
		if (del > 255) del = 255;
	}
	if (org->z < 55) {
		del = (5 * org->z + 245);
		if (del > 255) del = 255;
	}
	unsigned int Cr = GetColor(155, 55, 55);
	pr.z -= 20;
	pf.z += 20;
	p1.x += L * cos(PI / 180 * (*objang + org->zm));
	p1.y += L * sin(PI / 180 * (*objang + org->zm));
	p2.x += L * cos(PI / 180 * (*objang + 90 + org->zm));
	p2.y += L * sin(PI / 180 * (*objang + 90 + org->zm));
	p3.x += L * cos(PI / 180 * (*objang + 180 + org->zm));
	p3.y += L * sin(PI / 180 * (*objang + 180 + org->zm));
	p4.x += L * cos(PI / 180 * (*objang + 270 + org->zm));
	p4.y += L * sin(PI / 180 * (*objang + 270 + org->zm));
	if (org->active == 1) {
		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 220 - del);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y, org->hand[1], FALSE);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y, org->hand[1], FALSE);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y, org->hand[1], FALSE);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y, org->hand[1], FALSE);

		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 155 - del);
	}

	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p1, &vp).x, prs(&p1, &vp).y, Cr);
	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p2, &vp).x, prs(&p2, &vp).y, Cr);
	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p3, &vp).x, prs(&p3, &vp).y, Cr);
	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p4, &vp).x, prs(&p4, &vp).y, Cr);

	SetDrawBlendMode(DX_BLENDMODE_ALPHA, 255 - del);
	DrawRect(
		prs(&p1, &vp).x, prs(&p1, &vp).y,
		prs(&p2, &vp).x, prs(&p2, &vp).y,
		prs(&p3, &vp).x, prs(&p3, &vp).y,
		prs(&p4, &vp).x, prs(&p4, &vp).y,
		Cr, 1);

	if (org->active == 1) {
		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 100 - del);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y, org->hand[1], FALSE);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y, org->hand[1], FALSE);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y, org->hand[1], FALSE);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y, org->hand[1], FALSE);
		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 255 - del);
	}

	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p1, &vp).x, prs(&p1, &vp).y, Cr);
	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p2, &vp).x, prs(&p2, &vp).y, Cr);
	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p3, &vp).x, prs(&p3, &vp).y, Cr);
	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p4, &vp).x, prs(&p4, &vp).y, Cr);
	SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 0);

}
void DrawObject1(obj *org, double *objang) {
	obj p1 = *org, p2 = *org, p3 = *org, p4 = *org, pr = *org, pf = *org;
	int L = 100;
	double del = 0;
	if (org->z > 480) {
		del = 3 * (org->z - 480);
		if (del > 255) del = 255;
	}
	if (org->z < 55) {
		del = (5 * org->z + 245);
		if (del > 255) del = 255;
	}
	unsigned int Cr = GetColor(55, 55, 55);
	pr.z -= 20;
	pf.z += 20;
	p1.x += L * cos(PI / 180 * (*objang + org->zm));
	p1.y += L * sin(PI / 180 * (*objang + org->zm));
	p2.x += L * cos(PI / 180 * (*objang + 90 + org->zm));
	p2.y += L * sin(PI / 180 * (*objang + 90 + org->zm));
	p3.x += L * cos(PI / 180 * (*objang + 180 + org->zm));
	p3.y += L * sin(PI / 180 * (*objang + 180 + org->zm));
	p4.x += L * cos(PI / 180 * (*objang + 270 + org->zm));
	p4.y += L * sin(PI / 180 * (*objang + 270 + org->zm));
	if (org->active == 1) {
		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 220 - del);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y, org->hand[2], FALSE);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y, org->hand[2], FALSE);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y, org->hand[2], FALSE);
		DrawModiGraph(
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&pr, &vp).x, prs(&pr, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y, org->hand[2], FALSE);

		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 155 - del);
	}

	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p1, &vp).x, prs(&p1, &vp).y, Cr);
	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p2, &vp).x, prs(&p2, &vp).y, Cr);
	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p3, &vp).x, prs(&p3, &vp).y, Cr);
	DrawLine(
		prs(&pr, &vp).x, prs(&pr, &vp).y,
		prs(&p4, &vp).x, prs(&p4, &vp).y, Cr);

	SetDrawBlendMode(DX_BLENDMODE_ALPHA, 255 - del);
	DrawRect(
		prs(&p1, &vp).x, prs(&p1, &vp).y,
		prs(&p2, &vp).x, prs(&p2, &vp).y,
		prs(&p3, &vp).x, prs(&p3, &vp).y,
		prs(&p4, &vp).x, prs(&p4, &vp).y,
		Cr, 1);

	if (org->active == 1) {
		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 100 - del);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y, org->hand[2], FALSE);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p2, &vp).x, prs(&p2, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y, org->hand[2], FALSE);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p3, &vp).x, prs(&p3, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y, org->hand[2], FALSE);
		DrawModiGraph(
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&pf, &vp).x, prs(&pf, &vp).y,
			prs(&p4, &vp).x, prs(&p4, &vp).y,
			prs(&p1, &vp).x, prs(&p1, &vp).y, org->hand[2], FALSE);
		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 255 - del);
	}

	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p1, &vp).x, prs(&p1, &vp).y, Cr);
	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p2, &vp).x, prs(&p2, &vp).y, Cr);
	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p3, &vp).x, prs(&p3, &vp).y, Cr);
	DrawLine(
		prs(&pf, &vp).x, prs(&pf, &vp).y,
		prs(&p4, &vp).x, prs(&p4, &vp).y, Cr);
	SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 0);

}
void DrawObject2(obj *org, double *objang) {
#define G 6
#define R 12

	obj p[R][G];
	int L = 100, i, k;
	int a = org->z;

	double del = 0;
	if (org->z > 480) {
		del = 3 * (org->z - 480);
		if (del > 255) del = 255;
	}
	if (org->z < 55) {
		del = (5 * org->z + 245);
		if (del > 255) del = 255;
	}
	if (org->active == 0)L += (org->z - 405) * 5;
	unsigned int Cr = GetColor(0, 55, 255);
	for (i = 0; i < R; i++) {
		for (k = 0; k < G; k++) {
			p[i][k] = *org;
			p[i][k].z -= i * 20 / (R - 1) + 10;

			p[i][k].x += (L - i) * cos(PI / 180 * (org->z + 360 / G * k + i * 3));
			p[i][k].y += (L - i) * sin(PI / 180 * (org->z + 360 / G * k + i * 3));

		}
	}
	if (org->active == 1) {
	}
	SetDrawBlendMode(DX_BLENDMODE_ADD, 255 - del);
	for (i = 0; i < R; i++) {
		for (k = 0; k < G; k++) {
			DrawLine(
				prs(&p[i][k], &vp).x, prs(&p[i][k], &vp).y,
				prs(&p[i][(k + 1) % G], &vp).x, prs(&p[i][(k + 1) % G], &vp).y, Cr, 2);
		}
	}
	SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 0);

}
void DrawObject3(obj *org, double *objang) {
#define G 5
#define R 12

	obj p[R][G];
	int L = 100, i, k;
	int a = org->z;

	double del = 0;
	if (org->z > 480) {
		del = 3 * (org->z - 480);
		if (del > 255) del = 255;
	}
	if (org->z < 55) {
		del = (5 * org->z + 245);
		if (del > 255) del = 255;
	}
	unsigned int Cr = GetColor(55, 35, 1);
	if (org->active == 0)L += (org->z - 405)*0.5;

	for (i = 0; i < R; i++) {
		for (k = 0; k < G; k++) {
			p[i][k] = *org;
			p[i][k].z -= i * 20 / (R - 1) + 10;
			if (org->active == 0) {
				p[i][k].x += (L - i) * cos(PI / 180 * (-org->z*1.5 + 360 / G * k - i * 3));
				p[i][k].y += (L - i) * sin(PI / 180 * (-org->z*1.5 + 360 / G * k - i * 3));
			}
			if (org->active == 1) {
				p[i][k].x += (L - i) * cos(PI / 180 * (-org->z*0.3 + 360 / G * k - i * 3));
				p[i][k].y += (L - i) * sin(PI / 180 * (-org->z*0.3 + 360 / G * k - i * 3));
			}
		}
	}
	if (org->active == 1) {
	}
	SetDrawBlendMode(DX_BLENDMODE_ADD, 255 - del);
	for (i = 0; i < R; i++) {
		for (k = 0; k < G; k++) {
			DrawLine(
				prs(&p[i][k], &vp).x, prs(&p[i][k], &vp).y,
				prs(&p[i][(k + 1) % G], &vp).x, prs(&p[i][(k + 1) % G], &vp).y, Cr, 2);
		}
	}
	SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 0);

}
void DrawObject(obj *org, double *objang) {
	org->pat %= 4;
	if (org->pat == 0) {
		DrawObject0(org, objang);
	}
	if (org->pat == 1) {
		DrawObject1(org, objang);
	}
	if (org->pat == 2) {
		DrawObject2(org, objang);
	}
	if (org->pat == 3) {
		DrawObject3(org, objang);
	}
}

//ショットを描くよ
void DrawShot(obj *org, double *objang) {
	if (org->active == 1) {
		obj p[4][3];
		int i, k;

		double del = 0;
		if (org->z > 480) {
			del = 3 * (org->z - 480);
			if (del > 255) del = 255;
		}
		if (org->z < 250) {
			del = 2 * (328 - org->z);
			if (del > 255) del = 255;
			if (del < 0) del = 0;
		}
		for (i = 0; i < 4; i++) {
			p[i][2] = *org;
			p[i][2].x += 35 * cos(PI / 180 * i * 90);
			p[i][2].y += 35 * sin(PI / 180 * i * 90);
		}
		for (i = 0; i < 4; i++) {
			p[i][1] = *org;
			p[i][1].x += 40 * cos(PI / 180 * i * 90);
			p[i][1].z += 40 * sin(PI / 180 * i * 90);
		}
		for (i = 0; i < 4; i++) {
			p[i][0] = *org;
			p[i][0].z += 40 * cos(PI / 180 * i * 90);
			p[i][0].y += 40 * sin(PI / 180 * i * 90);
		}
		if (org->active == 1) {
		}
		SetDrawBlendMode(DX_BLENDMODE_ALPHA, 255 - del);
		for (k = 0; k < 3; k++) {
			DrawModiGraph(
				prs(&p[0][k], &vp).x, prs(&p[0][k], &vp).y,
				prs(&p[1][k], &vp).x, prs(&p[1][k], &vp).y,
				prs(&p[2][k], &vp).x, prs(&p[2][k], &vp).y,
				prs(&p[3][k], &vp).x, prs(&p[3][k], &vp).y, org->hand[0], TRUE);
		}
		SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 0);
		if (org->z < 200) {
			org->active = 0;
		}

	}
}

//距離感を計るためのよくわからない光る三角形を描くよ
void DrawSh(obj *org, double *ang) {
	obj p1 = *org, p2 = *org, p3 = *org;
	unsigned int Cr = GetColor(100, 155, 55);

	SetDrawBlendMode(DX_BLENDMODE_ADD, org->z*0.5);
	p1.x += 50 * cos(PI / 180 * (*ang + 0 + org->zm));
	p1.y += 50 * sin(PI / 180 * (*ang + 0 + org->zm));
	p2.x += 50 * cos(PI / 180 * (*ang + 120 + org->zm));
	p2.y += 50 * sin(PI / 180 * (*ang + 120 + org->zm));
	p3.x += 50 * cos(PI / 180 * (*ang + 240 + org->zm));
	p3.y += 50 * sin(PI / 180 * (*ang + 240 + org->zm));
	DrawTriangle(
		prs(&p1, &vp).x, prs(&p1, &vp).y,
		prs(&p2, &vp).x, prs(&p2, &vp).y,
		prs(&p3, &vp).x, prs(&p3, &vp).y, Cr, 0);
	SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);
}

//バイタリティトライアングルの表示
void tri(vpt org, unsigned int Cr, int hand) {
	vpt p[3];
	int i;
	for (i = 0; i < 3; i++) {
		p[i] = org;
		p[i].x += 20 * sin(PI / 180 * (counter + i * 120));
		p[i].y += 20 * cos(PI / 180 * (counter + i * 120));
	}

	for (i = 0; i < 3; i++) {
		DrawModiGraph(
			org.x, org.y,
			p[i].x, p[i].y,
			p[(i + 1) % 3].x, p[(i + 1) % 3].y,
			org.x, org.y, hand, FALSE);
		DrawLine(p[i].x, p[i].y, p[(i + 1) % 3].x, p[(i + 1) % 3].y, Cr, 2);

	}
	SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);

}

//インターフェースを描くよ
void Interface(int vital, int traceV, int bomb, double boost, double V, int UI[]) {

	vpt org;//orgはUI全体の座標の基準位置
	org.x = 0;
	org.y = 0;
	vpt vit[10];
	vpt bom = org;
	vpt bst;

	int i, k;
	int time;
	unsigned int Cr = GetColor(100, 155, 55);

	if (boost > 0) {
		if (boost >= 180) SetDrawBlendMode(DX_BLENDMODE_ADD, (1 + sin(PI / 180 * (450 - boost) * 5) * 40));
		if (boost<180) SetDrawBlendMode(DX_BLENDMODE_ADD, (1 + sin(PI / 180 * (450 - boost) * 10) * 40));
		DrawBox(0, 0, 640, 480, GetColor(
			255,
			255,
			255), 1);
		bst.x = prs(&pl, &vp).x;
		bst.y = prs(&pl, &vp).y;
		bst.y += 50;
		if (boost < 180) {
			if (boost>160)time = boost - 160;

			SetDrawBlendMode(DX_BLENDMODE_ALPHA, 240-12*time);
			DrawBox(bst.x - 30 * int(boost) / 180, bst.y, bst.x + 30 * int(boost) / 180, bst.y + 10, GetColor(
				255,
				0,
				0), 1);
			SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);
		}

	}
	//ブーストゲージの表示


	for (k = 0; k < traceV / 100; k++) {
		vit[k] = org;
		vit[k].x += 30 + 25 * k;
		vit[k].y += 50;
		SetDrawBlendMode(DX_BLENDMODE_ADD, 220 - Abs(100 * vital - traceV));
		tri(vit[k], Cr, UI[0]);
	}
	//体力ゲージの表示

	bom.x += 15;
	bom.y += 80;
	Cr = GetColor(0, 0, 0);
	DrawBox(bom.x, bom.y, bom.x + 130, bom.y + 25, Cr, 1);
	SetDrawBlendMode(DX_BLENDMODE_ADD, 255);
	Cr = GetColor((bomb - 10) / 7, 255 - (bomb - 10) / 7, 0);
	DrawBox(bom.x + 6, bom.y, bom.x + 6 + bomb / 15, bom.y + 25, Cr, 1);
	if (bomb == 1800) {
		SetDrawBlendMode(DX_BLENDMODE_ADD, 255);
		Cr = GetColor(
			255 * 0.5 * (1 + sin(PI / 180 * counter * 5)),
			255 * 0.5 * (1 + sin(PI / 180 * counter * 5)),
			255 * 0.5 * (1 + sin(PI / 180 * counter * 5)));
		DrawBox(bom.x, bom.y, bom.x + 130, bom.y + 25, Cr, 1);
	}
	//ボムゲージの表示
	SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);	DrawGraph(bom.x, bom.y, UI[1], 1);
}

struct DrawSet {
	void(*Draw)(obj*, double*);
	obj *org;
	double *ang;
};
#define Q 120 //障害物の数だよ
DrawSet srs[9 + Q];

bool Pre(DrawSet o1, DrawSet o2) {
	return o1.org->z < o2.org->z;
}

void DrawAll() {
	std::sort(srs, srs + 8 + Q, Pre);
	for (int i = 0; i<Q + 9; i++) {
		//		GraphFilterBlt(*srs[i].org->hand,DX_GRAPH_FILTER_GAUSS, 106, 100);
		if (srs[i].org->z < 700 && srs[i].org->z > 0)srs[i].Draw(srs[i].org, srs[i].ang);
		//z座標700以上0未満だと表示されなくなるよ
	}
}



struct objpat {
	int row;//xy平面一枚あたりに並ぶオブジェの数
	int column;//z座標の数
	int inter;//xy平面のz座標の間隔
	int rota;//回転の有無
	int rand;//0 or 1 位置の無作為選定 EV[].randと掛け合わせて使う
	int per[100] = { 0 };//それぞれのオブジェが出現する割合
	int req;//必要なオブジェクトの数(行×列)
};
struct objpat ord[50];
//パターンに応じて行列数間隔回転の有無を決めて構造体にまとめて返す
objpat Order(int pat) {
	int i=0, k=0;
	struct objpat a;
	int per[4] = { 0 }, cnt;

	if (pat == 0) {
		i = GetRand(5);
		if (i == 5)a.rota = 1;
		else a.rota = 0;
		a.row = GetRand(1) + 1;
		a.column = 30;
		a.inter = 150;
		a.rand = 1;
		per[0] = 85; per[1] = 10; per[2] = 4; per[3] = 1;
	}
	if (pat == 1) {
		i = GetRand(5);
		if (i == 5)a.row = 3;
		if (5 > i > 0)a.row = 2;
		if (i == 0)a.row = 1;

		a.column = 25;
		a.inter = 80;
		a.rand = 0;
		a.rota = GetRand(1);
		per[0] = 100-100*GetRand(1); per[1] = 100-per[0]; per[2] = 0; per[3] = 0;
	}
	if (pat == 2) {
		a.rota = 1 - GetRand(2);
		a.row = 1;
		a.column = 100;
		a.inter = 100;
		a.rand = 1;
		per[0] = 97; per[1] = 0; per[2] = 3; per[3] = 0;

	}
	if (pat == 3) {
		a.rota = 1 - GetRand(2);
		a.row = 2;
		a.column = 1;
		a.inter = 10;
		a.rand = 1;
		per[0] = 0; per[1] = 0; per[2] = 50; per[3] = 50;

	}
	if (pat == 4) {
		a.rota = 1 - GetRand(2);
		a.row = 5;
		a.column = 10;
		a.inter = 200;
		a.rand = 0;
		per[0] = 80; per[1] = 15; per[2] = 4; per[3] = 1;

	}
	if (pat == 5) {
		a.rota = 1 - GetRand(2);
		a.row = 3;
		a.column = 33;
		a.inter = 200;
		a.rand = 1;
		per[0] = 40; per[1] = 40; per[2] = 20; per[3] = 0;

	}
	if (pat == 6) {
		a.rota = 1;
		a.row = 3;
		a.column = 20;
		a.inter = 200;
		a.rand = 0;
		per[0] = 50; per[1] = 50; per[2] = 0; per[3] = 0;

	}
	a.req = a.row*a.column;
	i = 0;
	for (k = 0; k < 4; k++) {
		for (cnt = 0; cnt < per[k]; cnt++) {
			a.per[i] = k;
			i++;
		}
	}
	return a;
}

//ここからゲーム本体だよ
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {
	SetMainWindowText("Space Vector");
	ChangeWindowMode(1);
	SetFullScreenResolutionMode(DX_FSRESOLUTIONMODE_DESKTOP);
	DxLib_Init();
	SetDrawScreen(DX_SCREEN_BACK);

	FILE *fp;
	int scr = MakeScreen(640, 480, 0);
	int result;
	int i, k, j;
	int gamephase = 0;//ゲーム進行度
	int hittime = 0;//無敵時間
	int interval = 0;//ショットインターバル
	int ina = 20; //インターバル計算用
	int pat[Q] = { 0 };//出現パターン
	int vital = 5;//体力
	int traceV = vital * 100;//体力を追いかける
	int bomb = 1800;//ボムゲージ
	int traceB = 0;//ボムゲージを追いかける
	int flag;//障害物が出尽くしたかを判断
	int now, next;//現在の障害物の出現パターン
	double Erand[Q] = { 0 };
	double last;//障害物一群の最後尾
	double bombtime = 0;//ボム時間
	double vec = 0, vabs = 0;//移動関連
	double boost[2] = { 0,0 };//ブースト時間　加速分
	double ang = 90, objang = 0, shotang = ang;
	double V, VS = 3;
	double W;//自機が中心からどのくらい離れてるか(固定)
	unsigned int Cr;//色
	obj p[N], q[N], r[N], s[N], t[N];//こっちは描画順は気にしない
	obj  sh[10], EV[Q], shot[3];//こっちの変数は気にする
	void(*Draw)(obj*, double*);

	pl.z = 425;
	pl.pat = 0;

	//	
	for (i = 0; i < Q; i++){
		EV[i].active = 0;
	}
	now = 0;
	ord[0] = Order(0);
	i = 0;
	for (k = 0; k < ord[0].column; k++) {
		for (j = 0; j < ord[0].row; j++) {
			while (EV[i].active == 1) i++; 
			EV[i].z = -k * ord[0].inter;
			EV[i].active = 1;
			EV[i].pat = ord[0].per[GetRand(99)];
			Erand[i] = GetRand(359);
			i++;
		}
	}


	for (i = 0; i < N; i++) {
		p[i].z = i * 60;
		q[i].z = i * 60;
		r[i].z = i * 60;
		s[i].z = i * 60;
	}

	for (i = 0; i < Q; i++) {
		LoadDivGraph("pl.png", 3, 3, 1, 50, 50, EV[i].hand, 0);
	}
	int Music[10] = {
		LoadSoundMem("jingle.wav") ,
		LoadSoundMem("loop.mp3") ,
		LoadSoundMem("bom17_b.wav") ,
		LoadSoundMem("tm2_laser002.wav") ,
		LoadSoundMem("tm2_power001.wav") ,
		LoadSoundMem("tm2_shoot003.wav") ,
		LoadSoundMem("") ,
		LoadSoundMem("tm2_bom003.wav") ,
		LoadSoundMem("tm2_death000.wav") ,
		LoadSoundMem("tm2_bom002.wav") };
	int iq = LoadGraph("backobj.png");
	int title = LoadGraph("titlelogo.png");
	int pnrm = LoadGraph("pn.png");
	int road = LoadGraph("road.png");
	int UI[10] = { 0 };
	UI[0] = LoadGraph("vital.png");
	UI[1] = LoadGraph("bom.png");
	p[0].hand[0] = LoadGraph("player1.png");
	shot[0].hand[0] = LoadGraph("shot.png");
	shot[1].hand[0] = LoadGraph("shot.png");
	shot[2].hand[0] = LoadGraph("shot.png");
	LoadDivGraph("pl.png", 3, 3, 1, 50, 50, pl.hand, 0);
	pl.hand[3] = LoadGraph("boost.png");
	//画像ハンドル
#define M 4.5
	PlaySoundMem(Music[0], DX_PLAYTYPE_BACK, TRUE);

	while (!ProcessMessage()) {//処理開始
		ClearDrawScreen();
		if (gamephase == 0) {
			DrawGraph(0, 0, pnrm, FALSE);
			for (i = 0; i < N; i++) {
				p[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 135 + M * i));
				p[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 135 + M * i));
				s[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 45 + M * i));
				s[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 45 + M * i));
				q[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 225 + M * i));
				q[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 225 + M * i));
				r[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 315 + M * i));
				r[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 315 + M * i));
				p[i].z += 3; q[i].z += 3; r[i].z += 3; s[i].z += 3;
				if (p[i].z >= 1200) {
					p[i].z = 0; q[i].z = 0; r[i].z = 0; s[i].z = 0;
				}
			}
			for (i = 0; i < N; i++) {
				SetDrawBlendMode(DX_BLENDMODE_ALPHA, p[i].z / 2);
				DrawRect(prs(&p[i], &vp).x, prs(&p[i], &vp).y,
					prs(&q[i], &vp).x, prs(&q[i], &vp).y,
					prs(&r[i], &vp).x, prs(&r[i], &vp).y,
					prs(&s[i], &vp).x, prs(&s[i], &vp).y, GetColor(25, 235, 255), 1);
			}
			SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 1);
			DrawRotaGraph(320, 140, 1.5, 0, title, 1);
			DrawFormatString(260, 300, GetColor(255, 255, 255), "PRESS Z KEY");
			fopen_s(&fp, "highscore.dat", "rb");
			if (fp != NULL) {
				fread(&result, sizeof(result), 1, fp);
				fclose(fp);
				DrawFormatString(260, 250, GetColor(255, 255, 255), "HISCORE %d", result);
			}
			else {
				result = 0;
				fopen_s(&fp, "highscore.dat", "wb");
				if (fp != NULL) {
					fwrite(&result, sizeof(result), 1, fp);
					fclose(fp);
				}
				fopen_s(&fp, "highscore.dat", "rb");
				if (fp != NULL) {
					fread(&result, sizeof(result), 1, fp);
					fclose(fp);
					DrawFormatString(260, 250, GetColor(255, 255, 255), "HISCORE %d", result);
				}
			}

			ScreenFlip();
			if (CheckHitKey(KEY_INPUT_Z)) {
				PlaySoundMem(Music[1], DX_PLAYTYPE_LOOP, TRUE);
				gamephase = 1;
			}
		}
		if (gamephase == 1) {
			counter++;
			pl.zm++;
			//各座標のセット
			W = 150 + 1 * sin(PI / 12 * pl.zm);
			pl.x = 320 + 1.3* W * cos(PI / 180 * ang);
			pl.y = 240 + W * sin(PI / 180 * ang);
			vp.x = 320 + (pl.x - vp.x)*0.5; vp.y = 240 + (pl.y - vp.y)*0.5;
			V = VS + boost[1];
			//消失点の計算

			for (i = 0; i < 5; i++) {
				sh[i] = pl;
				sh[i].z -= i * 60 - 1;
			}
			//自機の距離感(ry
			for (i = 0; i < Q; i++) {
				EV[i].v = 0;
			}

			objang++;
			i = 0;
			for (k = 0; k < ord[now].column; k++) {
				for (j = 0; j < ord[now].row; j++) {
					if (EV[i].active == 1) {
						EV[i].x = 320 + 1.3*W * cos(PI / 180 *
							(Erand[i] * ord[now].rand + ord[now].rota * objang + 10 * k + j * 360 / ord[now].row));
						EV[i].y = 240 + W * sin(PI / 180 *
							(Erand[i] * ord[now].rand + ord[now].rota * objang + 10 * k + j * 360 / ord[now].row));
					}
					if (EV[i].z > 700) {
						EV[i].active = 0;
						EV[i].z = 700;
					}
					if ((EV[i].active == 0 && EV[i].z <= 0) == FALSE) {
						EV[i].z += V;
					}
					i++;
				}
			}
			flag = 0;
			for (i = 0; i < Q; i++) {
				if (700 > EV[i].z && EV[i].z > 0) {
					flag = 1;
				}
			}
			if (flag == 0) {
				Score += 10000;
				ord[now] = Order(GetRand(6));
				i = 0;
				for (k = 0; k < ord[0].column; k++) {
					for (j = 0; j < ord[0].row; j++) {
						while (EV[i].active == 1) i++;
						EV[i].z = -k * ord[0].inter;
						EV[i].active = 1;
						EV[i].pat = ord[0].per[GetRand(99)];
						Erand[i] = GetRand(359);
						i++;
					}
				}

			}
			//ここにつぎのパターン

			//オブジェクトの出るパターンとか

			for (i = 0; i < Q; i++) {
				EV[i].zm += 4;
				srs[i].Draw = DrawObject;
				srs[i].org = &EV[i];
				srs[i].ang = &objang;
			}
			for (i = Q; i < Q + 5; i++) {
				sh[i].zm += 4;
				srs[i].Draw = DrawSh;
				srs[i].org = &sh[i - Q];
				srs[i].ang = &ang;
			}
			srs[Q + 5].Draw = DrawPlayer;
			srs[Q + 5].org = &pl;
			srs[Q + 5].ang = &ang;
			for (i = 0; i < 3; i++) {
				srs[Q + 6 + i].Draw = DrawShot;
				srs[Q + 6 + i].org = &shot[i];
				srs[Q + 6 + i].ang = &shotang;
			}
			//ここで描画順を気にする奴らを関数アドレスに突っ込む


			for (i = 0; i < N; i++) {
				p[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 135 + M * i));
				p[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 135 + M * i));
				s[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 45 + M * i));
				s[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 45 + M * i));
				q[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 225 + M * i));
				q[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 225 + M * i));
				r[i].x = 320 + 400 * cos(PI / 180 * (ang * 0.5 + 315 + M * i));
				r[i].y = 240 + 400 * sin(PI / 180 * (ang * 0.5 + 315 + M * i));
				p[i].z += V; q[i].z += V; r[i].z += V; s[i].z += V;
				if (p[i].z >= 1200) {
					p[i].z = 0; q[i].z = 0; r[i].z = 0; s[i].z = 0;
				}
			}
			//↑背景ワイヤーのセット

			//入力関連
			ina = 20;
			if (vital > 0) {
				if (CheckHitKey(KEY_INPUT_Z)) {
					for (i = 0; i < 3; i++) {
						if (shot[i].active == 0 && interval == 0) {
							shot[i] = pl;
							shot[i].z -= 2.5 * 25;
							shot[i].hand[0] = LoadGraph("shot.png");
							shot[i].active = 1;
							interval += ina;
							break;
						}
					}
				}
				if (bombtime <= 0) {
					if (bomb >= 1800) {
						if (CheckHitKey(KEY_INPUT_X)) {
							PlaySoundMem(Music[7], DX_PLAYTYPE_BACK, TRUE);

							bomb = 0;
							bombtime = 60;
						}
					}
					else {
						bomb++;
					}
				}
				if (interval > 0) {
					if (interval == ina)pl.z += 5;
					if (interval == ina - 2)pl.z += 5;
					if (interval == ina - 4)pl.z += 5;
					if (interval == ina - 6)pl.z -= 3;
					if (interval == ina - 8)pl.z -= 3;
					if (interval == ina - 10)pl.z -= 3;
					if (interval == ina - 12)pl.z -= 3;
					if (interval == ina - 14)pl.z -= 3;
					interval--;
				}
				if (CheckHitKey(KEY_INPUT_RIGHT)) {
					vec -= 0.5;
					if (vec < -3.0)vec = -3.0;
				}
				if (CheckHitKey(KEY_INPUT_LEFT)) {
					vec += 0.5;
					if (vec > 3.0)vec = 3.0;
				}
				if (vec > 0)vec -= 0.125;
				if (vec < 0)vec += 0.125;
				ang += vec;
				if (ang < 0)ang += 360;
			}
			//ここから描画
			DrawGraph(0, 0, pnrm, FALSE);
			//遠景

			//ここになんかキラキラしたやつ

			for (i = 0; i < N; i++) {
				SetDrawBlendMode(DX_BLENDMODE_ALPHA, p[i].z / 2);
				DrawRect(prs(&p[i], &vp).x, prs(&p[i], &vp).y,
					prs(&q[i], &vp).x, prs(&q[i], &vp).y,
					prs(&r[i], &vp).x, prs(&r[i], &vp).y,
					prs(&s[i], &vp).x, prs(&s[i], &vp).y, GetColor(25, 235, 255), 1);
			}
			SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);
			//ペルソナ４のテレビ入るときのみたいなやつ

			//			objang+=2;
			for (k = 0; k < 3; k++) {
				if (shot[k].active == 1) {
					shot[k].z -= 7;
					for (i = 0; i < Q; i++) {
						if (EV[i].z - 10 < shot[k].z && shot[k].z < EV[i].z + 10
							&& (shot[k].x - EV[i].x)*(shot[k].x - EV[i].x) + (shot[k].y - EV[i].y)*(shot[k].y - EV[i].y) <= 110 * 110
							&& EV[i].active == 1
							&& EV[i].pat == 0) {
							Score += 500;
							EV[i].active = 0;
							shot[k].active = 0;
							//ショット命中
						}
					}
				}
			}
			if (hittime <= 0) {
				for (i = 0; i < Q; i++) {
					if (EV[i].z - 10 < pl.z && pl.z < EV[i].z + 10
						&& (pl.x - EV[i].x)*(pl.x - EV[i].x) + (pl.y - EV[i].y)*(pl.y - EV[i].y) <= 110 * 110
						&& EV[i].active == 1) {
						if (EV[i].pat == 0 || EV[i].pat == 1) {
							hittime = 10;
							if (pl.pat != 2) {
								PlaySoundMem(Music[2], DX_PLAYTYPE_BACK, TRUE);
								vital -= 1;

							}
							else {
								hittime += 10;
							}
							if (pl.pat == 2)Score += 500;
							//HP減
						}
						if (EV[i].pat == 2) {
							PlaySoundMem(Music[3], DX_PLAYTYPE_BACK, TRUE);
							boost[0] = 450;
							hittime += 10;
							pl.pat = 2;
							//加速モードオン
						}
						if (EV[i].pat == 3) {
							PlaySoundMem(Music[4], DX_PLAYTYPE_BACK, TRUE);
							vital = 5;
							bomb = 1800;
							//HPBOMB回復
						}
						EV[i].active = 0;
						//消滅時演出void関数をここに入れる animation(EV[i])

					}

				}
			}
			else {
				if (hittime > 0) {
					SetDrawBlendMode(DX_BLENDMODE_ALPHA, 20 * hittime);
					if (pl.pat != 2) {
						DrawBox(0, 0, 640, 480,
							GetColor(255, 0, 0), 1);
						vp.x += GetRand(50) - 25;
						vp.y += GetRand(50) - 25;
						pl.x += GetRand(50) - 25;
						pl.y += GetRand(50) - 25;

					}
					if (pl.pat == 2) {
						DrawBox(0, 0, 640, 480,
							GetColor(200, 200, 255), 1);
					}
					SetDrawBlendMode(DX_BLENDMODE_ALPHA, 255);
				}
				hittime--;
			}

			if (boost[0] > 0) {
				if (boost[0] > 290)boost[1] += 0.3;
				if (boost[1] > 3)boost[1] = 3;
				if (boost[0] <= 50)boost[1] -= 0.06;
				if (boost[0] == 1) {
					boost[1] = 0;
					pl.pat = 0;
				}
				boost[0] -= 1;
			}

			if (bombtime > 0) {
				for (i = 0; i < Q; i++) {
					if (EV[i].active == 1
						&& 600>EV[i].z
						&& EV[i].z > 50
						&& (EV[i].pat == 0 || EV[i].pat == 1)) {
						EV[i].active = 0;
						Score += 500;
					}
				}

				SetDrawBlendMode(DX_BLENDMODE_ADD, 4.25*bombtime);
				DrawRotaGraph(prs(&pl, &vp).x, prs(&pl, &vp).y, 121 - 2 * bombtime, 0, pl.hand[3], TRUE);
				vp.x += GetRand(2 * bombtime) - bombtime;
				vp.y += GetRand(2 * bombtime) - bombtime;
				bombtime--;
			}


			//当たり判定の計算は自機とオブジェクトの中心座標の距離で行う。（円形）


			DrawAll();
			//関数アドレスに突っ込まれたものは
			//ここでソートされてまとめて表示される

			if (bomb > traceB) {
				for (i = 0; i < 180; i++) {
					traceB += 1;
					if (bomb == traceB) {
						break;
					}
				}
			}
			if (bomb < traceB) {
				for (i = 0; i < 360; i++) {
					traceB -= 1;
					if (bomb == traceB) {
						break;
					}
				}
			}
			if (100 * vital > traceV) {
				traceV += (100 * vital - traceV) / 20 + 10;
			}
			if (100 * vital < traceV) {
				for (i = 0; i < 5; i++) {
					traceV -= 1;
					if (100 * vital == traceV) {
						break;
					}
				}
			}
			if (vital > 0) {
				Score += V;
			}
			if (Score > 50000 && VS < 4)VS += 1;
			if (Score > 100000 && VS < 5)VS += 1;
			if (Score > 150000 && VS < 6)VS += 1;
			if (Score > 200000 && VS < 7)VS += 1;
			if (Score > 250000 && VS < 8)VS += 1;
			if (Score > 300000 && VS < 9)VS += 1;
			DrawFormatString(10, 10, GetColor(255, 255, 255), "SCORE %10.0f", Score);
			Interface(vital, traceV, traceB, boost[0], V, UI);
			DrawFormatString(40, 84, GetColor(255, 255, 255), "BOMB %3d%%", bomb / 18);

			//スコアとかのインターフェース

			Cr = GetColor(255, 0, 0);

			if (vital <= 0) {
				if (vital == 0) {
					SetDrawBlendMode(DX_BLENDMODE_INVDESTCOLOR, 255);
					DrawBox(0, 0, 640, 480, GetColor(255, 255, 255), 1);
					ScreenFlip();
					StopSoundMem(Music[1]);
					StopSoundMem(Music[2]);
					PlaySoundMem(Music[8], DX_PLAYTYPE_BACK, TRUE);
					WaitTimer(800);
					PlaySoundMem(Music[9], DX_PLAYTYPE_BACK, TRUE);
				}
				SetDrawBlendMode(DX_BLENDMODE_ALPHA, 15-vital*4);
				DrawBox(prs(&pl, &vp).x+vital*12, 0, prs(&pl, &vp).x - vital * 12,480, GetColor(255, 255, 255),1);
				SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);
				DrawRotaGraph(prs(&pl, &vp).x, prs(&pl, &vp).y, -2 * vital + 1, 0, pl.hand[3], TRUE);
				DrawRotaGraph(prs(&pl, &vp).x, prs(&pl, &vp).y, -2 * vital + 1, 0, pl.hand[3], TRUE);
				vp.x += GetRand(2 * (60 + vital)) - (60 + vital);
				vp.y += GetRand(2 * (60 + vital)) - (60 + vital);
				vital--;
				SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);
				if (vital <= -60) {
					WaitTimer(1000);
					i = 0;
					fopen_s(&fp, "highscore.dat", "rb");
					if (fp != NULL) {
						fread(&result, sizeof(result), 1, fp);
						fclose(fp);
						if (int(Score) > result) {
							i = 1;
							result = int(Score);
							fopen_s(&fp, "highscore.dat", "wb");
							if (fp != NULL) {
								fwrite(&result, sizeof(result), 1, fp);
								fclose(fp);
							}
						}
					}
					gamephase = 2;
				}
			}
				
			if (counter < 15) {
				SetDrawBlendMode(DX_BLENDMODE_ALPHA, 255-counter*17);
				DrawBox(0, 0, 640, 480, 0, 1);
				SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 255);
			}
			ScreenFlip();

		}

		if (gamephase == 2) {
			k = 0;
			while (CheckHitKey(KEY_INPUT_ESCAPE)!=1){
				DrawBox(0, 0, 640, 480, GetColor(255-k,255-k,255-k), 1);
			if (i == 1) {
				DrawFormatString(220, 160, GetColor(k, 255-k, 255-k), "   NEW RECORD!!  ");
			}
			SetDrawBlendMode(DX_BLENDMODE_ALPHA, k);
			SetFontSize(22);
			DrawFormatString(220, 200, GetColor(k, k, k), "  SCORE   %7.0f", Score);
			DrawFormatString(220, 240, GetColor(k, k, k), "HIGHSCORE %7.0d", result);
			DrawFormatString(220, 280, GetColor(k, k, k), "  PRESS ESC KEY  ");
			ScreenFlip();
			k+=3;
			if(k>255)k=255;
			}
			break;
		}
	}

	DxLib_End();
	return 0;
}

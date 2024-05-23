#pragma once


class ViewArea {
public:
	ViewArea() {}
	ViewArea(double l, double b, double w, double h) : left(l), bottom(b), width(w), height(h) {}

	double left, bottom, width, height;
	const double aspectRatio() { return width / height; }
};
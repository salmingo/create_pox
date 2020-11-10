/**
 * create_pox 为ASA望远镜生成位置偏差文件
 * 输入:
 * - .fits：sequence生成的FITS文件
 * - .wcs: astrometry.net生成的WCS文件
 * 输出:
 * - PointError.pox
 **/

#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/param.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "FITSHandler.hpp"
#include "ADefine.h"
#include "ATimeSpace.h"

using namespace std;
using namespace AstroUtil;

ATimeSpace ats;
/*!
 * @struct one_point 单次天区指向位置统计信息
 */
struct one_point {
	string pathFITS;
	string pathWCS;
	string dateobs;
	string timeobs;
	double exptime;
	double objctra, objctdec;	// 当前历元
	int pierside;
	double crval1, crval2;		// J2000

public:
	/*!
	 * @brief 计算指向位置与天球位置的偏差
	 * @param era  赤经偏差, 量纲: 角度
	 * @param edc  赤纬偏差, 量纲: 角度
	 */
	void calc_error(double& era, double &edc) {
		if ((era = objctra - crval1) >= 12.0) era -= 24.0;
		else if (era < -12.0) era += 24.0;
		era *= 15.0;
		edc  = objctdec - crval2;
	}
};

void load_fits(one_point& pt) {
	FITSHandler hfit;
	if (hfit(pt.pathFITS.c_str())) {
		char dateobs[40], timeobs[40], *ptr;
		int notes;
		fits_read_key(hfit(), TSTRING, "DATE-OBS", dateobs,      NULL, &hfit.errcode);
		fits_read_key(hfit(), TDOUBLE, "EXPTIME",  &pt.exptime,  NULL, &hfit.errcode);
		fits_read_key(hfit(), TDOUBLE, "OBJCTRA",  &pt.objctra,  NULL, &hfit.errcode);
		fits_read_key(hfit(), TDOUBLE, "OBJCTDEC", &pt.objctdec, NULL, &hfit.errcode);
		fits_read_key(hfit(), TINT,    "NOTES",    &notes,       NULL, &hfit.errcode);
		if ((ptr = strstr(dateobs, "T")) != NULL)
			strcpy(timeobs, ptr + 1);

		pt.dateobs = dateobs;
		pt.timeobs = timeobs;
		pt.pierside = notes == 0 ? 1 : -1;
	}
}

void load_wcs(one_point& pt) {
	FITSHandler hfit;
	if (hfit(pt.pathWCS.c_str())) {
		fits_read_key(hfit(), TDOUBLE, "CRVAL1",  &pt.crval1,  NULL, &hfit.errcode);
		fits_read_key(hfit(), TDOUBLE, "CRVAL2",  &pt.crval2,  NULL, &hfit.errcode);
		pt.crval1 /= 15.0;
	}
}

void change_object(one_point& pt) {
	FITSHandler hfit;
	if (hfit(pt.pathFITS.c_str()), 1) {
		fits_update_key(hfit(), TDOUBLE, "OBJCTRA",  &pt.crval1, NULL, &hfit.errcode);
		fits_update_key(hfit(), TDOUBLE, "OBJCTDEC", &pt.crval2, NULL, &hfit.errcode);
	}
}

/*!
 * @brief 坐标从J2000转换至当前历元
 */
void epoch2current(one_point& pt) {
	int YY, MM, DD, hh, mm;
	double ss;
	double ra, dec;

	sscanf (pt.dateobs.c_str(), "%d-%d-%dT%d:%d:%lf",
			&YY, &MM, &DD,
			&hh, &mm, &ss);
	ats.SetUTC(YY, MM, DD, (hh + (mm + ss / 60.0) / 60.0) / 24.0);
	ats.EqTransfer(pt.crval1 * 15.0 * D2R, pt.crval2 * D2R, ra, dec);
	pt.crval1 = ra * R2H;
	pt.crval2 = dec * R2D;
}

int main(int argc, char **argv) {
	ats.SetSite(117.57444, 40.39583, 900.0, 8);

	char cwd[MAXPATHLEN];
	struct dirent **entry_list;
	struct dirent *entry;
	int count, n;
	vector<one_point> all_point;
	char pathwcs[40];
	char *ptr;

	// 统计当前目录下对应的FITS文件和WCS文件, 以及WCS文件的数量
	getcwd(cwd, MAXPATHLEN);
	count = scandir(cwd, &entry_list, 0, alphasort);
	for (int i = 0; i < count; ++i) {
		entry = entry_list[i];
		if ((ptr = strstr(entry->d_name, ".fit")) != NULL) {
			n = ptr - entry->d_name;
			strncpy(pathwcs, entry->d_name, n);
			strcpy(pathwcs + n, ".wcs");
			if (!access(pathwcs, R_OK)) {
				one_point pt;
				pt.pathFITS = entry->d_name;
				pt.pathWCS  = pathwcs;
				all_point.push_back(pt);
			}
		}

		free(entry);
	}
	free(entry_list);
	printf ("%d files found\n\n", count = all_point.size());

	// 生成目标文件
	double era, edc;
	double sum_ra(0.0), sum_dec(0.0), sq_ra(0.0), sq_dec(0.0);
	FILE *fp = fopen(argc == 1 ? "PointError.pox" : argv[1], "w");
	fprintf (fp, "%d\r\n", count);
	for (n = 0; n < count; ++n) {
		one_point& pt = all_point[n];
		// 从FITS文件提取DATE-OBS/EXPTIME/OBJCTRA/OBJCTDEC/NOTES
		// 从WCS文件提取CRVAL1和CRVAL2
		load_fits(pt);
		load_wcs(pt);
//		epoch2current(pt);	// (CRVAL1, CRVAL2)由J2000转换至当前历元
		pt.calc_error(era, edc);
		sum_ra  += era;
		sum_dec += edc;
		sq_ra   += era * era;
		sq_dec  += edc * edc;

		fprintf (fp,
				"\"Number %d\"\r\n"
				"\"\'%s\'\"\r\n"
				"\"%s\"\r\n"
				"\"%.3f\"\r\n"
				"%.6f\r\n"
				"%.6f\r\n"
				"%.6f\r\n"
				"%.6f\r\n"
				"\"%d\"\r\n"
				"\"*****************************\"\r\n",
				n + 1,
				pt.dateobs.c_str(),
				pt.timeobs.c_str(),
				pt.exptime,
				pt.objctra,
				pt.crval1,
				pt.objctdec,
				pt.crval2,
				pt.pierside);

		printf ("%s\n", pt.pathFITS.c_str());
		printf ("Number %d\n"
				"\tDATE-OBS = %s\n"
				"\tTIME-OBS = %s\n"
				"\tEXPTIME  = %.3f\n"
				"\tOBJCTRA  = %.5f\n"
				"\tCRVAL1   = %.5f\n"
				"\tOBJCTDEC = %.5f\n"
				"\tCRVAL2   = %.5f\n"
				"\tPIERSIDE = %d\n"
				"*****************************\n",
				n + 1,
				pt.dateobs.c_str(),
				pt.timeobs.c_str(),
				pt.exptime,
				pt.objctra,
				pt.crval1,
				pt.objctdec,
				pt.crval2,
				pt.pierside);
	}
	fclose(fp);

	printf ("\n\nPointing Error:\n");
	printf ("\tRA Error : %.3f arcsec\n"
			"\tDEC Error: %.3f arcsec\n\n",
			sqrt((sq_ra - sum_ra * sum_ra / count) / count) * 3600.0,
			sqrt((sq_dec - sum_dec * sum_dec / count) / count) * 3600.0
			);
	return 0;
}

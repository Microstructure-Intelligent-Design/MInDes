/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once
#include"time.h"
#include<sys/timeb.h>
#include<iostream>
#include<string>
namespace pf {
#define MESSURE
	using namespace std;
	struct ints_time {
		int year;
		int month;
		int day;
		int hour;
		int minute;
		int second;
		ints_time() {
			year = 0;
			month = 0;
			day = 0;
			hour = 0;
			minute = 0;
			second = 0;
		}
		ints_time& operator=(const ints_time& n) {
			year = n.year;
			month = n.month;
			day = n.day;
			hour = n.hour;
			minute = n.minute;
			second = n.second;
			return *this;
		}
		ints_time& operator+(const ints_time& n) {
			year += n.year;
			month += n.month;
			day += n.day;
			hour += n.hour;
			minute += n.minute;
			second += n.second;
			return *this;
		}
		ints_time& operator-(const ints_time& n) {
			year -= n.year;
			month -= n.month;
			day -= n.day;
			hour -= n.hour;
			minute -= n.minute;
			second -= n.second;
			return *this;
		}
		bool operator>(const ints_time& n) {
			if (year > n.year)
				return true;
			else if (year < n.year)
				return false;
			if (month > n.month)
				return true;
			else if (month < n.month)
				return false;
			if (day > n.day)
				return true;
			else if (day < n.day)
				return false;
			if (hour > n.hour)
				return true;
			else if (hour < n.hour)
				return false;
			if (minute > n.minute)
				return true;
			else if (minute < n.minute)
				return false;
			if (second > n.second)
				return true;
			else if (second < n.second)
				return false;
			return false;
		}
		bool operator==(const ints_time& n) {
			if (year == n.year && month == n.month && day == n.day && hour == n.hour && minute == n.minute && second == n.second)
				return true;
			return false;
		}
	};
	namespace timer{
		static void init(double& total_time_begin);
		static double total_duration_sec(double& total_time_begin);
		static double total_duration_min(double& total_time_begin);
		static double total_duration_hour(double& total_time_begin);
		static void print_cunrrent_time_on_screen(); 
		static string return_cunrrent_time_by_string();
		static string get_cunrrent_time_as_file_name();
		static void get_cunrrent_time_by_int(int& year, int& month, int& day, int& hour, int& minute, int& second);
		static void get_cunrrent_time_by_ints_time(ints_time& t);
		static void interval_begin(double& interval_begin);
		static double interval_end(double& interval_begin);  //By Second
		static void time_interval_precision_secs_begin(double& t_interval);
		static void time_interval_precision_secs_end(double& t_interval);
		static string trans_int_time_to_string(ints_time time);
	};
	inline double timer::interval_end(double& interval_begin) {
		timeb t;
		ftime(&t);//获取毫秒
		double _now = t.time + t.millitm * 0.001;
		return _now - interval_begin;
	}
	inline void timer::interval_begin(double& interval_begin) {
		timeb t;
		ftime(&t);//获取毫秒
		interval_begin = t.time + t.millitm * 0.001;
	}
	inline void timer::init(double& total_time_begin){
		timeb t;
		ftime(&t);//获取毫秒
		total_time_begin = t.time + t.millitm * 0.001;
	}
	inline double timer::total_duration_sec(double& total_time_begin){
		timeb t;
		ftime(&t);//获取毫秒
		double _now = t.time + t.millitm * 0.001;
		return _now - total_time_begin;
	}
	inline double timer::total_duration_min(double& total_time_begin) {
		timeb t;
		ftime(&t);//获取毫秒
		double _now = t.time + t.millitm * 0.001;
		return (_now - total_time_begin) / 60.0;
	}
	inline double timer::total_duration_hour(double& total_time_begin) {
		timeb t;
		ftime(&t);//获取毫秒
		double _now = t.time + t.millitm * 0.001;
		return (_now - total_time_begin) / 3600.0;
	}
	inline string timer::trans_int_time_to_string(ints_time time) {
		return to_string(time.year) + "-" + to_string(time.month) + "-" + to_string(time.day)
			+ "-" + to_string(time.hour) + "-" + to_string(time.minute) + "-" + to_string(time.second);
	}
	inline void timer::print_cunrrent_time_on_screen() {
		time_t timer;
		time(&timer);
#ifdef _WIN32
		char buf[26];
		ctime_s(buf, 26, &timer);                                     ///< Windows style directory separator
#else
		char* buf;
		buf = ctime(&timer);                                         ///< Unix/Linux style directory separator
#endif
		std::cout << "## " << buf;
	}
	inline string timer::return_cunrrent_time_by_string() {
		time_t timer;
		time(&timer);
#ifdef _WIN32
		char buf[26];
		ctime_s(buf, 26, &timer);                                     ///< Windows style directory separator
#else
		char* buf;
		buf = ctime(&timer);                                         ///< Unix/Linux style directory separator
#endif
		string time(buf);
		return "## " + time;
	}
	inline string timer::get_cunrrent_time_as_file_name() {
		time_t timer;
		time(&timer);
#ifdef _WIN32
		char buf[26];
		ctime_s(buf, 26, &timer);                                     ///< Windows style directory separator
#else
		char* buf;
		buf = ctime(&timer);                                         ///< Unix/Linux style directory separator
#endif
		string time(buf);

		for (auto charr = time.begin(); charr < time.end();) {
			if (*charr == '\n')
				time.erase(charr);
			else if (*charr == ':') {
				*charr = '-';
				charr++;
			}
			else
				charr++;
		}

		return time;
	}
	inline void timer::get_cunrrent_time_by_int(int& year, int& month, int& day, int& hour, int& minute, int& second) {
		const time_t t = time(NULL);
		struct tm* system_time = localtime(&t);
		year = 1900 + system_time->tm_year;
		month = 1 + system_time->tm_mon;
		day = system_time->tm_mday;
		hour = system_time->tm_hour;
		minute = system_time->tm_min;
		second = system_time->tm_sec;
	}
	inline void timer::get_cunrrent_time_by_ints_time(ints_time& tt) {
		const time_t t = time(NULL);
		struct tm* system_time = localtime(&t);
		tt.year = 1900 + system_time->tm_year;
		tt.month = 1 + system_time->tm_mon;
		tt.day = system_time->tm_mday;
		tt.hour = system_time->tm_hour;
		tt.minute = system_time->tm_min;
		tt.second = system_time->tm_sec;
	}
	static void timer::time_interval_precision_secs_begin(double& t_interval) {
		timeb t;
		ftime(&t);//获取毫秒
		t_interval = t.time + t.millitm * 0.001;
	}
	static void timer::time_interval_precision_secs_end(double& t_interval) {
		timeb t;
		ftime(&t);//获取毫秒
		t_interval = (t.time + t.millitm * 0.001) - t_interval;
	}
}
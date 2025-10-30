#pragma once
#include <chrono>
#include <iostream>
#include <string>
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
	namespace timer {
		inline void init(double& total_time_begin) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			total_time_begin = ticks * seconds_per_tick;
			// auto time_point_now = std::chrono::system_clock::now();
			// total_time_begin = time_point_now.time_since_epoch().count() / REAL(1e9);
		}
		inline double total_duration_sec(double& total_time_begin) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			double _now = ticks * seconds_per_tick;
			return _now - total_time_begin;
		}
		inline double total_duration_min(double& total_time_begin) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			double _now = ticks * seconds_per_tick;
			return (_now - total_time_begin) / double(60.0);
		}
		inline double total_duration_hour(double& total_time_begin) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			double _now = ticks * seconds_per_tick;
			return (_now - total_time_begin) / double(3600.0);
		}
		inline double time_min(double& sec_time) {
			return sec_time / double(60.0);
		}
		inline double time_hour(double& sec_time) {
			return sec_time / double(3600.0);
		}
		inline void print_cunrrent_time_on_screen() {
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
		inline string return_cunrrent_time_by_string() {
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
		inline string get_cunrrent_time_as_file_name() {
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
		inline void get_cunrrent_time_by_int(int& year, int& month, int& day, int& hour, int& minute, int& second) {
			const time_t t = time(NULL);
			struct tm system_time;
#ifdef _WIN32
			localtime_s(&system_time, &t);
#else
			localtime_r(&t, &system_time);
#endif
			year = 1900 + system_time.tm_year;
			month = 1 + system_time.tm_mon;
			day = system_time.tm_mday;
			hour = system_time.tm_hour;
			minute = system_time.tm_min;
			second = system_time.tm_sec;
		}
		inline void get_cunrrent_time_by_ints_time(ints_time& tt) {
			const time_t t = time(NULL);
			struct tm system_time;
#ifdef _WIN32
			localtime_s(&system_time, &t);
#else
			localtime_r(&t, &system_time);
#endif
			tt.year = 1900 + system_time.tm_year;
			tt.month = 1 + system_time.tm_mon;
			tt.day = system_time.tm_mday;
			tt.hour = system_time.tm_hour;
			tt.minute = system_time.tm_min;
			tt.second = system_time.tm_sec;
		}
		inline void interval_begin(double& interval_begin) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			interval_begin = ticks * seconds_per_tick;
			//auto time_point_now = std::chrono::system_clock::now();
			//interval_begin = REAL(time_point_now.time_since_epoch().count() / 1e9);
		}
		inline double interval_end(double& interval_begin) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			double _now = ticks * seconds_per_tick;
			//auto time_point_now = std::chrono::system_clock::now();
			//REAL _now = REAL(time_point_now.time_since_epoch().count() / 1e9);
			return _now - interval_begin;
		}  //By Second
		inline void time_interval_precision_secs_begin(double& t_interval) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			t_interval = ticks * seconds_per_tick;
			//auto time_point_now = std::chrono::system_clock::now();
			//t_interval = REAL(time_point_now.time_since_epoch().count() / 1e9);
		}
		inline void time_interval_precision_secs_end(double& t_interval) {
			auto ticks = std::chrono::system_clock::now().time_since_epoch().count();
			double seconds_per_tick = (double)std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;
			t_interval = ticks * seconds_per_tick - t_interval;
			//auto time_point_now = std::chrono::system_clock::now();
			//t_interval = REAL(time_point_now.time_since_epoch().count() / 1e9) - t_interval;
		}
		inline string trans_int_time_to_string(ints_time time) {
			return to_string(time.year) + "-" + to_string(time.month) + "-" + to_string(time.day)
				+ "-" + to_string(time.hour) + "-" + to_string(time.minute) + "-" + to_string(time.second);
		}
	};
}
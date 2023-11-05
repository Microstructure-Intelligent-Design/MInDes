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
#include "../Base.h"

namespace pf {
	namespace CpuMemoryUsage {
        // get current process pid
        static int GetCurrentPid() {
            //return _getpid();
#ifdef WIN32
            return _getpid();
#else
            return getpid();
#endif
        }
        // get specific process cpu occupation ratio by pid
#ifdef WIN32
        static uint64_t convert_time_format(const FILETIME* ftime) {
            LARGE_INTEGER li;
            li.LowPart = ftime->dwLowDateTime;
            li.HighPart = ftime->dwHighDateTime;
            return li.QuadPart;
        }
#else
// FIXME: can also get cpu and mem status from popen cmd
// the info line num in /proc/{pid}/status file
#define VMRSS_LINE 22
#define PROCESS_ITEM 14
        static const char* get_items(const char* buffer, unsigned int item) {
            // read from buffer by offset
            const char* p = buffer;
            int len = strlen(buffer);
            int count = 0;
            for (int i = 0; i < len; i++)
            {
                if (' ' == *p)
                {
                    count++;
                    if (count == item - 1)
                    {
                        p++;
                        break;
                    }
                }
                p++;
            }
            return p;
        }
        static unsigned long get_cpu_total_occupy() {
            // get total cpu use time
            // different mode cpu occupy time
            unsigned long user_time;
            unsigned long nice_time;
            unsigned long system_time;
            unsigned long idle_time;
            FILE* fd;
            char buff[1024] = { 0 };
            fd = fopen("/proc/stat", "r");
            if (nullptr == fd)
                return 0;
            fgets(buff, sizeof(buff), fd);
            char name[64] = { 0 };
            sscanf(buff, "%s %ld %ld %ld %ld", name, &user_time, &nice_time, &system_time, &idle_time);
            fclose(fd);
            return (user_time + nice_time + system_time + idle_time);
        }
        static unsigned long get_cpu_proc_occupy(int pid) {
            // get specific pid cpu use time
            unsigned int tmp_pid;
            unsigned long utime;  // user time
            unsigned long stime;  // kernel time
            unsigned long cutime; // all user time
            unsigned long cstime; // all dead time
            char file_name[64] = { 0 };
            FILE* fd;
            char line_buff[1024] = { 0 };
            sprintf(file_name, "/proc/%d/stat", pid);
            fd = fopen(file_name, "r");
            if (nullptr == fd)
                return 0;
            fgets(line_buff, sizeof(line_buff), fd);
            sscanf(line_buff, "%u", &tmp_pid);
            const char* q = get_items(line_buff, PROCESS_ITEM);
            sscanf(q, "%ld %ld %ld %ld", &utime, &stime, &cutime, &cstime);
            fclose(fd);
            return (utime + stime + cutime + cstime);
        }
#endif
        static float GetCpuUsageRatio(int pid) {
#ifdef WIN32
            static int64_t last_time = 0;
            static int64_t last_system_time = 0;
            FILETIME now;
            FILETIME creation_time;
            FILETIME exit_time;
            FILETIME kernel_time;
            FILETIME user_time;
            int64_t system_time;
            int64_t time;
            int64_t system_time_delta;
            int64_t time_delta;
            // get cpu num
            SYSTEM_INFO info;
            GetSystemInfo(&info);
            int cpu_num = info.dwNumberOfProcessors;
            int cpu_ratio = 0;
            // get process hanlde by pid
            HANDLE process = OpenProcess(PROCESS_ALL_ACCESS, FALSE, pid);
            // use GetCurrentProcess() can get current process and no need to close handle
            // get now time
            GetSystemTimeAsFileTime(&now);
            if (!GetProcessTimes(process, &creation_time, &exit_time, &kernel_time, &user_time)) {
                // We don't assert here because in some cases (such as in the Task Manager)  
                // we may call this function on a process that has just exited but we have  
                // not yet received the notification.  
                printf("GetCpuUsageRatio GetProcessTimes failed\n");
                return 0.0;
            }
            // should handle the multiple cpu num
            system_time = (convert_time_format(&kernel_time) + convert_time_format(&user_time)) / cpu_num;
            time = convert_time_format(&now);
            if ((last_system_time == 0) || (last_time == 0)) {
                // First call, just set the last values.  
                last_system_time = system_time;
                last_time = time;
                return 0.0;
            }
            system_time_delta = system_time - last_system_time;
            time_delta = time - last_time;
            CloseHandle(process);
            if (time_delta == 0) {
                printf("GetCpuUsageRatio time_delta is 0, error\n");
                return 0.0;
            }
            // We add time_delta / 2 so the result is rounded.  
            cpu_ratio = (int)((system_time_delta * 100 + time_delta / 2) / time_delta); // the % unit
            last_system_time = system_time;
            last_time = time;
            return float(cpu_ratio / 100.0); // convert to float number
#else
            unsigned long totalcputime1, totalcputime2;
            unsigned long procputime1, procputime2;
            totalcputime1 = get_cpu_total_occupy();
            procputime1 = get_cpu_proc_occupy(pid);
            // FIXME: the 200ms is a magic number, works well
            usleep(200000); // sleep 200ms to fetch two time point cpu usage snapshots sample for later calculation
            totalcputime2 = get_cpu_total_occupy();
            procputime2 = get_cpu_proc_occupy(pid);
            float pcpu = 0.0;
            if (0 != totalcputime2 - totalcputime1)
                pcpu = (procputime2 - procputime1) / float(totalcputime2 - totalcputime1); // float number
            int cpu_num = get_nprocs();
            pcpu *= cpu_num; // should multiply cpu num in multiple cpu machine
            return pcpu;
#endif
        }
        // get specific process physical memeory occupation size by pid (MB)
        static double GetMemoryUsage(int pid) {
#ifdef WIN32
            uint64_t mem = 0, vmem = 0;
            PROCESS_MEMORY_COUNTERS pmc;
            // get process hanlde by pid
            HANDLE process = OpenProcess(PROCESS_ALL_ACCESS, FALSE, pid);
            if (GetProcessMemoryInfo(process, &pmc, sizeof(pmc))) {
                mem = pmc.WorkingSetSize;
                vmem = pmc.PagefileUsage;
            }
            CloseHandle(process);
            // use GetCurrentProcess() can get current process and no need to close handle
            // convert mem from B to MB
            return mem / 1024.0 / 1024.0;
#else
            char file_name[64] = { 0 };
            FILE* fd;
            char line_buff[512] = { 0 };
            sprintf(file_name, "/proc/%d/status", pid);
            fd = fopen(file_name, "r");
            if (nullptr == fd)
                return 0;
            char name[64];
            int vmrss = 0;
            for (int i = 0; i < VMRSS_LINE - 1; i++)
                fgets(line_buff, sizeof(line_buff), fd);
            fgets(line_buff, sizeof(line_buff), fd);
            sscanf(line_buff, "%s %d", name, &vmrss);
            fclose(fd);
            // cnvert VmRSS from KB to MB
            return vmrss / 1024.0;
#endif
        }
        static void init(FieldStorage_forPhaseNode& phaseMesh) {
            stringstream report;
            int current_pid = GetCurrentPid(); // or you can set a outside program pid
            //float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
            double memory_usage = GetMemoryUsage(current_pid);
            report << "> current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
            Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
        }
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
            stringstream report;
			int current_pid = GetCurrentPid(); // or you can set a outside program pid
			//float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
            double memory_usage = GetMemoryUsage(current_pid);
            report << "> current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
            Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
            stringstream report;
            int current_pid = GetCurrentPid(); // or you can set a outside program pid
            //float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
            double memory_usage = GetMemoryUsage(current_pid);
            report << "> current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
            return report.str();
		}
        static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
            stringstream report;
            int current_pid = GetCurrentPid(); // or you can set a outside program pid
            //float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
            double memory_usage = GetMemoryUsage(current_pid);
            report << "> current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
            Solvers::get_instance()->writer.add_string_to_txt_and_screen(report.str(), LOG_FILE_NAME);
        }
	}
}
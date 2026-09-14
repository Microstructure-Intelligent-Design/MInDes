#include "CpuMemoryUsage.h"
#include <cstdlib>
#include <cstring>

namespace pf {
    namespace cpu_memory_usage {
        // get current process pid
        int GetCurrentPid() {
            //return _getpid();
#ifdef _WIN32
            return _getpid();
#else
            return getpid();
#endif
        }
        // get specific process cpu occupation ratio by pid
// FIXME: can also get cpu and mem status from popen cmd
// the info line num in /proc/{pid}/status file
#define VMRSS_LINE 22
#define PROCESS_ITEM 14
#ifdef _WIN32
        uint64_t convert_time_format(const FILETIME* ftime) {
            LARGE_INTEGER li;
            li.LowPart = ftime->dwLowDateTime;
            li.HighPart = ftime->dwHighDateTime;
            return li.QuadPart;
        }
#else
        const char* get_items(const char* buffer, unsigned int item) {
            // read from buffer by offset
            const char* p = buffer;
            const std::size_t len = std::strlen(buffer);
            int count = 0;
            for (std::size_t i = 0; i < len; i++)
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
        unsigned long get_cpu_total_occupy() {
            // get total cpu use time
            // different mode cpu occupy time
            unsigned long user_time = 0;
            unsigned long nice_time = 0;
            unsigned long system_time = 0;
            unsigned long idle_time = 0;
            FILE* fd;
            char buff[1024] = { 0 };
            fd = fopen("/proc/stat", "r");
            if (nullptr == fd)
                return 0;
            if (fgets(buff, sizeof(buff), fd) == nullptr) {
                fclose(fd);
                return 0;
            }
            char name[64] = { 0 };
            if (sscanf(buff, "%63s %lu %lu %lu %lu", name, &user_time, &nice_time, &system_time, &idle_time) != 5) {
                fclose(fd);
                return 0;
            }
            fclose(fd);
            return (user_time + nice_time + system_time + idle_time);
        }
        unsigned long get_cpu_proc_occupy(int pid) {
            // get specific pid cpu use time
            unsigned int tmp_pid = 0;
            unsigned long utime = 0;  // user time
            unsigned long stime = 0;  // kernel time
            unsigned long cutime = 0; // all user time
            unsigned long cstime = 0; // all dead time
            char file_name[64] = { 0 };
            FILE* fd;
            char line_buff[1024] = { 0 };
            sprintf(file_name, "/proc/%d/stat", pid);
            fd = fopen(file_name, "r");
            if (nullptr == fd)
                return 0;
            if (fgets(line_buff, sizeof(line_buff), fd) == nullptr) {
                fclose(fd);
                return 0;
            }
            if (sscanf(line_buff, "%u", &tmp_pid) != 1) {
                fclose(fd);
                return 0;
            }
            const char* q = get_items(line_buff, PROCESS_ITEM);
            if (sscanf(q, "%lu %lu %lu %lu", &utime, &stime, &cutime, &cstime) != 4) {
                fclose(fd);
                return 0;
            }
            fclose(fd);
            return (utime + stime + cutime + cstime);
        }
#endif
        float GetCpuUsageRatio(int pid) {
#ifdef _WIN32
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
            int cpu_num = sysconf(_SC_NPROCESSORS_ONLN);
            pcpu *= cpu_num; // should multiply cpu num in multiple cpu machine
            return pcpu;
#endif
        }
        double GetMemoryUsage(int pid) {
#ifdef _WIN32
            uint64_t mem = 0;
            PROCESS_MEMORY_COUNTERS pmc;

            // 
            HANDLE process = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, pid);
            if (process == nullptr) {
                return 0.0; // 
            }

            if (GetProcessMemoryInfo(process, &pmc, sizeof(pmc))) {
                mem = pmc.WorkingSetSize; // 
            }

            CloseHandle(process);
            return static_cast<double>(mem) / 1024.0 / 1024.0; // 

#else
            char file_name[256];
            FILE* fd;
            char line[512];
            int vmrss = 0;

            // 
            snprintf(file_name, sizeof(file_name), "/proc/%d/status", pid);
            fd = fopen(file_name, "r");
            if (!fd) {
                return 0.0; // 
            }

            // 
            while (fgets(line, sizeof(line), fd)) {
                if (std::strncmp(line, "VmRSS:", 6) == 0) {
                    // 
                    char* value_str = line + 6;
                    vmrss = std::atoi(value_str); // 
                    break;
                }
            }
            fclose(fd);

            return static_cast<double>(vmrss) / 1024.0; // 
#endif
        }
        void exec_pre_iii() {
            std::stringstream report;
            int current_pid = GetCurrentPid(); // or you can set a outside program pid
            //float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
            double memory_usage = GetMemoryUsage(current_pid);
            report << "# current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
            WriteLog(report.str());
        }
        void exec_post_ii() {
            if (show_loop_information::screen_output_step == 0)
                return;
            if (main_iterator::Current_ITE_step % show_loop_information::screen_output_step == 0) {
                std::stringstream report;
                int current_pid = GetCurrentPid(); // or you can set a outside program pid
                //float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
                double memory_usage = GetMemoryUsage(current_pid);
                report << "# current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
                WriteLog(report.str());
            }
        }
        void deinit() {
            std::stringstream report;
            int current_pid = GetCurrentPid(); // or you can set a outside program pid
            //float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
            double memory_usage = GetMemoryUsage(current_pid);
            report << "# current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
            WriteLog(report.str());
        }

        void init_cpu_memory_usage() {
            std::stringstream report;
            int current_pid = GetCurrentPid(); // or you can set a outside program pid
            //float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
            double memory_usage = GetMemoryUsage(current_pid);
            report << "# current memory usage: " << memory_usage << " MB ( " << memory_usage / 1024.0 << " GB )" << std::endl;
            WriteLog(report.str());
            load_a_new_module(nullptr, nullptr, exec_pre_iii,
                nullptr, nullptr, nullptr,
                nullptr, exec_post_ii, nullptr, deinit);
        }
    }
}
